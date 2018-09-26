package main

import (
  "bufio"
	"errors"
	"os"
	"time"

	"github.com/bmeg/grip/util"
	"github.com/cespare/xxhash"
	log "github.com/sirupsen/logrus"
	"github.com/spf13/cobra"
	"github.com/steakknife/bloomfilter"
)

var (
	maxElements uint64  = 100000000
	probCollide float64 = 0.00000001
	vertexFiles = []string{}
  vertexManifest = ""
	edgeFiles = []string{}
  edgeManifest = ""
	verbose   = false
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

    if vertexManifest != "" {
      verts, err := readLines(vertexManifest)
      if err != nil {
        return err
      }
      vertexFiles = append(vertexFiles, verts...)
    }
		if len(vertexFiles) == 0 {      
      return errors.New("must provide one or more vertex files to load")
    }

    if edgeManifest != "" {
      edges, err := readLines(edgeManifest)
      if err != nil {
        return err
      }
      edgeFiles = append(edgeFiles, edges...)
    }
		if len(edgeFiles) == 0 {
        return errors.New("must provide one or more edge files to check")
		}

		return check(vertexFiles, edgeFiles)
	},
}

func init() {
  log.SetFormatter(&log.JSONFormatter{TimestampFormat: time.RFC3339})
	flags := rootCmd.Flags()
	flags.StringSliceVarP(&vertexFiles, "vertex-file", "v", vertexFiles, "vertex file to load")
  flags.StringVarP(&vertexManifest, "vertex-manifest", "V", vertexManifest, "vertex file manifest")
	flags.StringSliceVarP(&edgeFiles, "edge-file", "e", edgeFiles, "edge file to check")
  flags.StringVarP(&edgeManifest, "edge-manifest", "E", edgeManifest, "edge file manifest")
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

func check(vertexFiles, edgeFiles []string) error {
	bf, err := bloomfilter.NewOptimal(maxElements, probCollide)
	if err != nil {
		return err
	}

	count := 0
	for _, f := range vertexFiles {
    log.Debugf("Loading vertices from: %s", f)
		for v := range util.StreamVerticesFromFile(f) {
			h := xxhash.New()
			h.Write([]byte(v.Gid))
			bf.Add(h)
			count++
			if count%1000 == 0 {
				log.Debugf("Loaded %d vertices", count)
			}
		}
	}
	log.Debugf("Loaded %d vertices", count)

	count = 0
	for _, f := range edgeFiles {
    log.Debugf("Loading edges from: %s", f)
		for e := range util.StreamEdgesFromFile(f) {
			f := xxhash.New()
			f.Write([]byte(e.From))
			if !bf.Contains(f) {
				log.WithFields(log.Fields{"Gid": e.Gid, "From": e.From}).Error("From references a vertex that does not exist")
			}
			t := xxhash.New()
			t.Write([]byte(e.To))
			if !bf.Contains(t) {
				log.WithFields(log.Fields{"Gid": e.Gid, "To": e.To}).Error("To references a vertex that does not exist")
			}
			count++
			if count%1000 == 0 {
				log.Debugf("Checked %d edges", count)
			}
		}
	}
	log.Debugf("Checked %d edges", count)

	return nil
}

func main() {
	if err := rootCmd.Execute(); err != nil {
		os.Exit(1)
	}
}
