package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"net/http"
	"os"
	"path/filepath"
)

type Interactions struct {
	Meta    Meta     `json:"_meta,omitempty"`
	Records []Record `json:"records,omitempty"`
}

type Meta struct {
	CurrentPage int32 `json:"current_page,omitempty"`
	// PerPage     int32 `json:"per_page,omitempty"`
	TotalPages int32 `json:"total_pages,omitempty"`
	Links      struct {
		Next     string `json:"next,omitempty"`
		Previous string `json:"previous,omitempty"`
	} `json:"links,omitempty"`
}

type Record struct {
	ID                string             `json:"id,omitempty"`
	GeneName          string             `json:"gene_name,omitempty"`
	EntrezID          int32              `json:"entrez_id,omitempty"`
	DrugName          string             `json:"drug_name,omitempty"`
	ChemblID          string             `json:"chembl_id,omitempty"`
	Publications      []int32            `json:"publications,omitempty"`
	InteractionTypes  []string           `json:"interaction_types,omitempty"`
	Sources           []string           `json:"sources,omitempty"`
	Attributes        []Attribute        `json:"attributes,omitempty"`
	InteractionClaims []InteractionClaim `json:"interaction_claims,omitempty"`
}

type Attribute struct {
	Name    string   `json:"name,omitempty"`
	Value   string   `json:"value,omitempty"`
	Sources []string `json:"sources,omitempty"`
}

type InteractionClaim struct {
	Source          string      `json:"source,omitempty"`
	Drug            string      `json:"drug,omitempty"`
	Gene            string      `json:"gene,omitempty"`
	IntractionTypes []string    `json:"interaction_types,omitempty"`
	Attributes      []Attribute `json:"attributes,omitempty"`
}

func main() {
	outputFile := ""
	flag.StringVar(&outputFile, "output", outputFile, "output file path")
	flag.Parse()

	var out io.WriteCloser
	var err error
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

	url := "http://dgidb.org/api/v2/interactions?count=100&page=1"
	for {
		resp, err := http.Get(url)
		if err != nil {
			panic(err)
		}
		defer resp.Body.Close()
		body, err := ioutil.ReadAll(resp.Body)
		if err != nil {
			panic(err)
		}

		interactions := Interactions{}
		err = json.Unmarshal(body, &interactions)
		if err != nil {
			fmt.Println(string(body))
			panic(err)
		}

		writer := json.NewEncoder(out)
		for _, i := range interactions.Records {
			err = writer.Encode(i)
			if err != nil {
				panic(err)
			}
		}

		if interactions.Meta.Links.Next == "" {
			break
		}
		url = interactions.Meta.Links.Next
	}
}
