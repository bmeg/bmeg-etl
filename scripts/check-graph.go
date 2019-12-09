package main

import (
	"bufio"
	"errors"
	"os"
	"strings"
	"time"

	"github.com/bmeg/grip/util"
	"github.com/cespare/xxhash"
	"github.com/gammazero/workerpool"
	log "github.com/sirupsen/logrus"
	"github.com/spf13/cobra"
	"github.com/steakknife/bloomfilter"
)

var (
	maxElements uint64  = 500000000
	probCollide float64 = 0.00000001
	manifest            = ""
	vertexFiles         = []string{}
	edgeFiles           = []string{}
	verbose             = false
	workers     int     = 10
)

var rootCmd = &cobra.Command{
	Use:   "check-graph",
	Short: "",
	Long:  ``,
	Args:  cobra.ExactArgs(0),
	RunE: func(cmd *cobra.Command, args []string) error {
		if verbose {
			log.SetLevel(log.DebugLevel)
		}

		files, err := readLines(manifest)
		if err != nil {
			return err
		}

		for _, fname := range files {
			if strings.HasSuffix(fname, ".Vertex.json.gz") {
				vertexFiles = append(vertexFiles, fname)
			} else if strings.HasSuffix(fname, ".Edge.json.gz") {
				edgeFiles = append(edgeFiles, fname)
			}
		}

		if len(vertexFiles) == 0 {
			return errors.New("must provide one or more vertex files in the manifest")
		}

		if len(edgeFiles) == 0 {
			return errors.New("must provide one or more edge files in the manifest")
		}

		return check(vertexFiles, edgeFiles, workers)
	},
}

func init() {
	log.SetFormatter(&log.JSONFormatter{TimestampFormat: time.RFC3339})
	flags := rootCmd.Flags()
	flags.StringVarP(&manifest, "manifest", "m", manifest, "file manifest listing vertex and edge files")
	flags.IntVarP(&workers, "workers", "w", workers, "number of workers to use to read the files and load bloom filter")
	flags.BoolVar(&verbose, "verbose", verbose, "verbose mode")
}

// readLines reads a whole file into memory
// and returns a slice of its lines.
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

func check(vertexFiles, edgeFiles []string, workers int) error {
	bf, err := bloomfilter.NewOptimal(maxElements, probCollide)
	if err != nil {
		return err
	}

	wp := workerpool.New(workers)
	for _, f := range vertexFiles {
		f := f
		wp.Submit(func() {
			log.WithFields(log.Fields{"file": f}).Debugf("Loading vertices")
			vChan, err := util.StreamVerticesFromFile(f)
			if err != nil {
				log.WithFields(log.Fields{"file": f}).Error(err)
				return
			}
			count := 0
			for v := range vChan {
				h := xxhash.New()
				_, err := h.Write([]byte(v.Gid))
				if err != nil {
					log.Error("xxhash.Write", err)
					continue
				}
				bf.Add(h)
				count++
				if count%1000 == 0 {
					log.WithFields(log.Fields{"file": f}).Debugf("Loaded %d vertices", count)
				}
			}
			log.WithFields(log.Fields{"file": f}).Debugf("Loaded %d vertices", count)
		})
	}
	wp.StopWait()

	wp = workerpool.New(workers)
	for _, f := range edgeFiles {
		f := f
		wp.Submit(func() {
			log.WithFields(log.Fields{"file": f}).Debugf("Loading edges")
			eChan, err := util.StreamEdgesFromFile(f)
			if err != nil {
				log.WithFields(log.Fields{"file": f}).Error(err)
				return
			}
			notFound := 0
			count := 0
			for e := range eChan {
				f := xxhash.New()
				_, err := f.Write([]byte(e.From))
				if err != nil {
					log.Error("xxhash.Write", err)
					continue
				}
				if !bf.Contains(f) {
					notFound++
					log.WithFields(log.Fields{"Gid": e.Gid, "From": e.From, "To": e.To, "Label": e.Label}).Error("From does not exist")
				}
				t := xxhash.New()
				_, err = t.Write([]byte(e.To))
				if err != nil {
					log.Error("xxhash.Write", err)
					continue
				}
				if !bf.Contains(t) {
					notFound++
					log.WithFields(log.Fields{"Gid": e.Gid, "From": e.From, "To": e.To, "Label": e.Label}).Error("To does not exist")
				}
				count++
				if count%1000 == 0 {
					log.Debugf("Checked %d edges", count)
				}
			}
			log.WithFields(log.Fields{"file": f}).Debugf("Checked %d edges", count)
			log.WithFields(log.Fields{"file": f}).Debugf("%d From / To references were not found", notFound)
		})
	}
	wp.StopWait()

	return nil
}

func main() {
	if err := rootCmd.Execute(); err != nil {
		os.Exit(1)
	}
}
