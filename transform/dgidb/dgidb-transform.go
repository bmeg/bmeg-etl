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
	// "github.com/golang/protobuf/jsonpb"
)

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

// CompoundIDs represents a subset of mappings from:
// https://www.ebi.ac.uk/unichem/rest/src_compound_id/{compound_id}/{source_id}
//
// Sources described here:
// https://www.ebi.ac.uk/unichem/ucquery/listSources
type CompoundID struct {
	// source_id 1
	ChEMBL string `json:"chembl,omitempty"`
	// source_id 22
	PubChem string `json:"pubchem,omitempty"`
	// source_id 2
	DrugBank string `json:"drugbank,omitempty"`
	// source_id 7
	ChEBI string `json:"chebi,omitempty"`
}

func main() {
	inputFile := ""
	outputFile := ""
	flag.StringVar(&inputFile, "input", inputFile, "input file")
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

	file, err := os.Open(inputFile)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	// writer := json.NewEncoder(out)
	// scanner := bufio.NewScanner(file)
	// for scanner.Scan() {
	// 	interaction := Record{}
	// 	err = json.Unmarshal(scanner.Bytes(), &interactions)
	// 	if err != nil {
	// 		panic(err)
	// 	}

	// 	err = writer.Encode(cid)
	// 	if err != nil {
	// 		panic(err)
	// 	}
	// }
}
