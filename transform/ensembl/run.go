package main

import (
	"bmeg"
	"compress/gzip"
	"flag"
	"fmt"
	"github.com/blachlylab/gff3"
	"github.com/golang/protobuf/jsonpb"
	"github.com/golang/protobuf/proto"
	"io"
	"log"
	"os"
	"strings"
)

type msg struct {
	name string
	msg  proto.Message
}

func id_clean(id string) string {

	for _, p := range []string{"exon:", "gene:", "transcript:"} {
		if strings.HasPrefix(id, p) {
			return id[len(p):len(id)]
		}
	}
	return id
}

func main() {
	flag.Parse()

	json := jsonpb.Marshaler{}

	out := make(chan msg)

	go func() {
		files := map[string]*os.File{}

		for m := range out {
			fname := fmt.Sprintf("%s.json", m.name)
			if _, ok := files[fname]; !ok {
				f, _ := os.Create(fname)
				files[fname] = f
			}
			s, _ := json.MarshalToString(m.msg)
			files[fname].Write([]byte(s))
			files[fname].Write([]byte("\n"))
		}

		for _, v := range files {
			v.Close()
		}
	}()

	gz_handle, _ := os.Open(flag.Arg(0))
	handle, _ := gzip.NewReader(gz_handle)
	read := gff3.NewReader(handle)
	for feat, err := read.Read(); err != io.EOF; feat, err = read.Read() {
		switch feat.TypeField {
		case "CDS":
			//fmt.Printf("%s %s\n", feat.AttributesField["ID"], feat.TypeField)
		case "five_prime_UTR":

		case "three_prime_UTR":
			//test

		case "biological_region":

		case "mRNA":
			t := bmeg.Transcript{
				Id:     id_clean(feat.AttributesField["ID"]),
				Parent: id_clean(feat.AttributesField["Parent"]),
				Symbol: feat.AttributesField["Name"],
				Strand: string(feat.StrandField),
				Start:  int32(feat.StartField),
				End:    int32(feat.EndField),
				SeqId:  feat.SeqidField,
			}
			out <- msg{"Transcript", &t}

		case "tRNA":
		case "vaultRNA_primary_transcript":
		case "lnc_RNA":

		case "gene":
			g := bmeg.Gene{
				Id:          id_clean(feat.AttributesField["ID"]),
				Symbol:      feat.AttributesField["Name"],
				Description: feat.AttributesField["description"],
				Strand:      string(feat.StrandField),
				Start:       int32(feat.StartField),
				End:         int32(feat.EndField),
				SeqId:       feat.SeqidField,
			}
			out <- msg{"Gene", &g}
		case "transcript":
			t := bmeg.Transcript{
				Id:     id_clean(feat.AttributesField["ID"]),
				Parent: id_clean(feat.AttributesField["Parent"]),
				Symbol: feat.AttributesField["Name"],
				Strand: string(feat.StrandField),
				Start:  int32(feat.StartField),
				End:    int32(feat.EndField),
				SeqId:  feat.SeqidField,
			}
			out <- msg{"Transcript", &t}
		case "processed_transcript":
		case "exon":
			t := bmeg.Exon{
				Id:     id_clean(feat.AttributesField["exon_id"]),
				Parent: id_clean(feat.AttributesField["Parent"]),
				Strand: string(feat.StrandField),
				Start:  int32(feat.StartField),
				End:    int32(feat.EndField),
				SeqId:  feat.SeqidField,
			}
			out <- msg{"Exon", &t}
		case "miRNA":
		case "snoRNA":
		case "lincRNA":
		case "lincRNA_gene":
		case "aberrant_processed_transcript":
		case "snRNA":
		case "rRNA":
		case "pseudogene":
		case "scRNA":
		case "chromosome":
		case "pseudogenic_transcript":
		case "processed_pseudogene":
		case "ncRNA":
		case "ncRNA_gene":
		case "V_gene_segment":
		case "J_gene_segment":
		case "C_gene_segment":
		case "D_gene_segment":
		case "VD_gene_segment":
		case "NMD_transcript_variant":
		case "snRNA_gene":
		case "miRNA_gene":
		case "RNA":
		case "snoRNA_gene":
		case "rRNA_gene":
		case "nc_primary_transcript":
		case "mt_gene":

		default:
			log.Printf("Unknown %s\n", feat.TypeField)
		}
	}
	gz_handle.Close()

	close(out)

}
