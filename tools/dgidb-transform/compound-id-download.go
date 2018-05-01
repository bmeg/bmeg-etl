package main

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"net/http"
	"os"
	"path/filepath"
)

// CompoundID represents a subset of mappings from:
// https://www.ebi.ac.uk/unichem/rest/src_compound_id/{compound_id}/{source_id}
//
// Sources described here:
// https://www.ebi.ac.uk/unichem/ucquery/listSources
func getCompoundIDs(chemblID string, srcMap map[string]string) (map[string]string, error) {
	compound := map[string]string{"chembl": chemblID}

	urlTmpl := "https://www.ebi.ac.uk/unichem/rest/src_compound_id/%s/1"
	respMap, err := httpGet(fmt.Sprintf(urlTmpl, chemblID))	
	if err != nil {
		return compound, err
	}

	for _, v := range respMap {
		compound[srcMap[v["src_id"]]] = v["src_compound_id"]
	}
	
	urlTmpl = "https://www.ebi.ac.uk/unichem/rest/structure/%s/1"
	respMap, err = httpGet(fmt.Sprintf(urlTmpl, chemblID))	
	if err != nil {
		return compound, err
	}

	if len(respMap) != 1 {
		return compound, fmt.Errorf("unexpected response from %s; %v", fmt.Sprintf(urlTmpl, chemblID), respMap) 
	}

	compound["standardinchi"] = respMap[0]["standardinchi"]
	compound["standardinchikey"] = respMap[0]["standardinchikey"]

	return compound, nil
}

func makeSourceMap() (map[string]string, error) {
	// src_id -> name
	srcMap := map[string]string{}

	srcURL := "https://www.ebi.ac.uk/unichem/rest/src_ids/"
	srcInfoURLTmpl := "https://www.ebi.ac.uk/unichem/rest/sources/%s"

	respMap, err := httpGet(srcURL)
	if err != nil {
		return nil, err
	}

	for _, src := range respMap {
		srcInfo, err := httpGet(fmt.Sprintf(srcInfoURLTmpl, src["src_id"]))	
		if err != nil {
			return nil, err
		}
		srcMap[src["src_id"]] = srcInfo[0]["name"]
	}

	return srcMap, nil
}

func httpGet(url string) ([]map[string]string, error) {
	resp, err := http.Get(url)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	body, err := ioutil.ReadAll(resp.Body)
	if err != nil {
		return nil, err
	}

	if resp.StatusCode != 200 {
		return nil, fmt.Errorf("[STATUS CODE - %d]\t%s", resp.StatusCode, body)
	}

	respMap := []map[string]string{}
	err = json.Unmarshal(body, &respMap)
	if err != nil {
		return nil, err
	}

	return respMap, nil
}

func main() {
	inputFile := ""
	outputFile := ""
	flag.StringVar(&inputFile, "input", inputFile, "input file containing a ChEMBL ID per line")
	flag.StringVar(&outputFile, "output", outputFile, "output file path")
	flag.Parse()

	if inputFile == "" {
		fmt.Println("interactions file must be provided")
		os.Exit(1)
	}

	file, err := os.Open(inputFile)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	var out io.WriteCloser
	if outputFile != "" {
		outputFile, err = filepath.Abs(outputFile)
		if err != nil {
			panic(err)
		}
		d := filepath.Dir(outputFile)
		err = os.MkdirAll(d, 0755)
		if err != nil {
			panic(err)
		}
		out, err = os.Create(outputFile)
		if err != nil {
			panic(err)
		}
	} else {
		out = os.Stdout
	}
	defer out.Close()

	srcMap, err := makeSourceMap()
	if err != nil {
		panic(err)
	}

	logger := log.New(os.Stderr, "logger: ", log.Lshortfile)

	writer := json.NewEncoder(out)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		chemblID := scanner.Text()
		cid, err := getCompoundIDs(chemblID, srcMap)
		if err != nil {
			logger.Print(err)
		}
		err = writer.Encode(cid)
		if err != nil {
			logger.Print(err)
		}
	}
}
