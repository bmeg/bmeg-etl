package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"github.com/biostream/schemas/go/bmeg"
	"github.com/golang/protobuf/jsonpb"
	"io"
	"os"
	"strconv"
	"strings"
)

type row struct {
	gene  string
	exprs []float64
}

func columns(s *bufio.Scanner) []string {
	return strings.Split(s.Text(), "\t")
}

func main() {
	filepath := ""
	source := ""
	scaleStr := "RPKM"
	gzipped := false

	flag.StringVar(&filepath, "filepath", filepath, "path to input file, in GCT format")
	flag.StringVar(&source, "source", source, "name of source, e.g. CCLE")
	flag.StringVar(&scaleStr, "scale", scaleStr, "expression scale, one of [RPKM, TPKM, READ_COUNT, FPKM]")
	flag.BoolVar(&gzipped, "gzipped", gzipped, "is the input file gzipped?")
	flag.Parse()

	if filepath == "" {
		fmt.Println("-filepath is required")
		os.Exit(1)
	}

	if source == "" {
		fmt.Println("-source is required")
		os.Exit(1)
	}

	if scaleStr == "" {
		fmt.Println("-scale is required")
		os.Exit(1)
	}

	msg := bmeg.GeneExpression{
		Source: source,
	}

	scaleInt, ok := bmeg.ExpressionScale_value[scaleStr]
	if !ok {
		panic("unrecognized scale")
	}
	msg.Scale = bmeg.ExpressionScale(scaleInt)

	f, err := os.Open(filepath)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	var reader io.Reader
	reader = f

	if gzipped {
		reader, err = gzip.NewReader(f)
		if err != nil {
			panic(err)
		}
	}

	scanner := bufio.NewScanner(reader)
	buf := make([]byte, 1024*512)
	scanner.Buffer(buf, 1024*512)

	// Drop version header
	scanner.Scan()

	// Load dimensions header
	scanner.Scan()
	dimensions := columns(scanner)
	rowCount, _ := strconv.Atoi(dimensions[0])
	columnCount, _ := strconv.Atoi(dimensions[1])

	// Load samples header
	var samples []string
	if scanner.Scan() {
		samples = columns(scanner)[2:]
	} else {
		panic(fmt.Errorf("no data found %s", scanner.Err()))
	}

	fmt.Fprintf(os.Stderr, "Loading data. genes %d, samples %d\n", rowCount, columnCount)

	// Load all the data into memory.
	var rows []row
	rowI := 0
	for scanner.Scan() {
		cols := columns(scanner)
		rowI++

		if rowI%1000 == 0 {
			fmt.Fprintf(os.Stderr, "Loaded %d rows of total %d\r", rowI, rowCount)
		}

		row := row{gene: cols[0]}
		for _, rawexpr := range cols[2:] {
			expr, err := strconv.ParseFloat(rawexpr, 64)
			if err != nil {
				panic(err)
			}
			row.exprs = append(row.exprs, expr)
		}
		rows = append(rows, row)
	}

	if err := scanner.Err(); err != nil {
		panic(err)
	}

	out, err := os.Create("out.json")
	if err != nil {
		panic(err)
	}
	defer out.Close()

	mar := jsonpb.Marshaler{}
	for sampleI, sample := range samples {
		fmt.Fprintf(os.Stderr, "Processing %s [%d of %d]\r", sample, sampleI, columnCount)
		msg.BiosampleId = sample
		msg.Expressions = make(map[string]float64, len(samples))

		for _, row := range rows {
			v := row.exprs[sampleI]
			if v == 0 {
				continue
			}
			msg.Expressions[row.gene] = v
		}
		mar.Marshal(out, &msg)
		out.Write([]byte("\n"))
	}
}
