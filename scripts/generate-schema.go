package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path"
	"regexp"
	"strings"
	"sync"
	"time"

	"github.com/bmeg/grip/gripql"
	"github.com/bmeg/grip/protoutil"
	"github.com/bmeg/grip/util"
	"github.com/gammazero/workerpool"
	"github.com/golang/protobuf/jsonpb"
	"github.com/iancoleman/strcase"
	log "github.com/sirupsen/logrus"
	"github.com/spf13/cobra"
)

var (
	manifest      = ""
	graphName     = ""
	sampleN   int = 100
	workers   int = 10
	yaml          = false
	verbose       = false
)

var rootCmd = &cobra.Command{
	Use:          "generate-schema",
	Short:        "",
	Long:         ``,
	SilenceUsage: true,
	Args:         cobra.ExactArgs(0),
	RunE: func(cmd *cobra.Command, args []string) error {
		if verbose {
			log.SetLevel(log.DebugLevel)
		} else {
			log.SetLevel(log.InfoLevel)
		}

		if manifest == "" {
			return fmt.Errorf("must provide a manifest file.\n\n%s", cmd.UsageString())
		}

		files, err := readLines(manifest)
		if err != nil {
			return err
		}

		if len(files) == 0 {
			return fmt.Errorf("must provide one or more vertex/edge files to load")
		}

		vChan := make(chan *vertex, 100)
		eChan := make(chan *edge, 100)
		wp := workerpool.New(workers)

		vertRE, err := regexp.Compile("(.*\\.)?(.*)(\\.Vertex.json.gz$)")
		if err != nil {
			return err
		}
		edgeRE, err := regexp.Compile("(.*\\.)?(.*)_(.*)_(.*)(\\.Edge.json.gz$)")
		if err != nil {
			return err
		}

		for _, fname := range files {
			fname := fname
			wp.Submit(func() {
				objSchema, err := getSchema(fname, sampleN)
				if err != nil {
					log.Errorf("error: %v", err)
					return
				}
				if strings.HasSuffix(fname, ".Vertex.json.gz") {
					match := vertRE.FindStringSubmatch(path.Base(fname))
					if len(match) != 4 {
						log.Error("vertex match error", path.Base(fname))
						return
					}
					label := match[2]
					vChan <- &vertex{label: label, data: objSchema}
				} else if strings.HasSuffix(fname, ".Edge.json.gz") {
					match := edgeRE.FindStringSubmatch(path.Base(fname))
					if len(match) != 6 {
						log.Error("edge match error", path.Base(fname))
						return
					}
					from := match[2]
					label := strcase.ToSnake(match[3])
					to := match[4]
					eChan <- &edge{
						label: label,
						from:  from,
						to:    to,
						data:  objSchema,
					}
				}
				return
			})
		}

		var waitgroup sync.WaitGroup
		waitgroup.Add(2)

		vertexMap := map[string]*vertex{}
		go func() {
			for v := range vChan {
				if vert, ok := vertexMap[v.label]; ok {
					data := util.MergeMaps(vert.data, v.data)
					v.data = data.(map[string]interface{})
				}
				vertexMap[v.label] = v
			}
			waitgroup.Done()
		}()

		edgeMap := map[edgeKey]*edge{}
		go func() {
			for e := range eChan {
				k := edgeKey{to: e.to, from: e.from, label: e.label}
				if edge, ok := edgeMap[k]; ok {
					data := util.MergeMaps(edge.data, e.data)
					e.data = data.(map[string]interface{})
				}
				edgeMap[k] = e
			}
			waitgroup.Done()
		}()

		// wait for all workers to finish
		log.Info("waiting for workers to complete")
		for {
			if wp.WaitingQueueSize() != 0 {
				log.Infof("queued tasks: %v", wp.WaitingQueueSize())
				time.Sleep(2 * time.Second)
			} else {
				log.Infof("queued tasks: %v", wp.WaitingQueueSize())
				break
			}
		}
		close(vChan)
		close(eChan)
		log.Info("workers finished")

		log.Info("waiting for schema processors to complete")
		waitgroup.Wait()
		log.Info("schema processors finished")

		vList := []*gripql.Vertex{}
		for label, v := range vertexMap {
			vSchema := &gripql.Vertex{
				Gid:   label,
				Label: label,
				Data:  protoutil.AsStruct(v.data),
			}
			vList = append(vList, vSchema)
		}

		eList := []*gripql.Edge{}
		for k, v := range edgeMap {
			eSchema := &gripql.Edge{
				Gid:   fmt.Sprintf("(%s)-%s->(%s)", k.from, k.label, k.to),
				Label: k.label,
				From:  k.from,
				To:    k.to,
				Data:  protoutil.AsStruct(v.data),
			}
			eList = append(eList, eSchema)
		}

		schema := &gripql.Graph{Graph: graphName, Vertices: vList, Edges: eList}

		var txt string
		if yaml {
			txt, err = gripql.GraphToYAMLString(schema)
		} else {
			txt, err = gripql.GraphToJSONString(schema)
		}
		if err != nil {
			return err
		}
		fmt.Printf("%s\n", txt)
		return nil
	},
}

func init() {
	log.SetFormatter(&log.JSONFormatter{TimestampFormat: time.RFC3339})
	flags := rootCmd.Flags()
	flags.StringVarP(&manifest, "manifest", "m", manifest, "file manifest listing vertex and edge files")
	flags.StringVarP(&graphName, "graph-name", "n", graphName, "name of the graph")
	flags.IntVarP(&sampleN, "sample", "s", sampleN, "number of elements to sample from each file")
	flags.BoolVar(&yaml, "yaml", yaml, "output schema in YAML rather than JSON format")
	flags.IntVarP(&workers, "workers", "w", workers, "number of workers to use to read the files and determine the schema")
	flags.BoolVarP(&verbose, "verbose", "v", verbose, "turn on verbose logging")
}

func main() {
	if err := rootCmd.Execute(); err != nil {
		os.Exit(1)
	}
}

func streamGzipLines(path string, n int) (chan []byte, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	r, err := gzip.NewReader(file)
	if err != nil {
		return nil, err
	}
	out := make(chan []byte, 100)
	go func() {
		defer file.Close()
		defer close(out)
		reader := bufio.NewReaderSize(r, 102400)
		var isPrefix bool = true
		var err error = nil
		var line, ln []byte
		i := 0
		for err == nil && i <= n {
			line, isPrefix, err = reader.ReadLine()
			ln = append(ln, line...)
			if !isPrefix {
				out <- ln
				ln = []byte{}
				i++
			}
		}
	}()
	return out, nil
}

func readLines(path string) ([]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var lines []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}
	return lines, scanner.Err()
}

func getSchema(file string, n int) (map[string]interface{}, error) {
	log.WithFields(log.Fields{"file": file}).Debug("processing")
	if strings.Contains(file, "Expression.Vertex.json") || strings.Contains(file, "CopyNumberAlteration.Vertex.json") {
		n = 1
		log.WithFields(log.Fields{"file": file}).Debugf("override n to %v", n)
	}
	reader, err := streamGzipLines(file, n)
	if err != nil {
		return nil, err
	}
	m := jsonpb.Unmarshaler{AllowUnknownFields: true}
	schema := map[string]interface{}{}
	i := 0
	for line := range reader {
		// all verts should be able to be unmarshalled into an Edge struct
		obj := &gripql.Edge{}
		err = m.Unmarshal(bytes.NewReader(line), obj)
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("unmarshaling: %v", err)
		}
		data := protoutil.AsMap(obj.Data)
		ds := gripql.GetDataFieldTypes(data)
		util.MergeMaps(schema, ds)
		i++
		if i%50 == 0 {
			log.WithFields(log.Fields{"file": file}).Debugf("processed %v / %v records", i, n)
		}
	}
	return schema, nil
}

type edgeKey struct {
	label, to, from string
}

type edge struct {
	label, to, from string
	data            map[string]interface{}
}

type vertex struct {
	label string
	data  map[string]interface{}
}
