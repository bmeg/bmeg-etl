package org.pcextract;

import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.io.BioPAXIOHandler;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.pattern.miner.*;

import java.io.*;
import java.net.URL;
import java.util.Set;
import java.util.zip.GZIPInputStream;

public class Main {
    public static void main(String[] args) throws IOException {

        String input = args[0];
        String outputs = args[1];

        SimpleIOHandler ioHandler = new SimpleIOHandler();
        FileInputStream fin = new FileInputStream(input);

        InputStream gin = new GZIPInputStream(fin);
        BioPAXIOHandler handler = new SimpleIOHandler();
        Model model = handler.convertFromOWL(gin);

        FileWriter outSIF = new FileWriter(outputs + ".extSIF");
        FileWriter outComplex = new FileWriter(outputs + ".complex");

        PrintWriter sifPW = new PrintWriter(outSIF);

        SIFSearcher sifSearch = new SIFSearcher(null, SIFEnum.values());
        Set<SIFInteraction> sifSet = sifSearch.searchSIF(model);
        for (SIFInteraction curLink : sifSet) {
            if (curLink.type != SIFEnum.NEIGHBOR_OF) {
                URL surl = new URL(curLink.sourceID);
                URL turl = new URL(curLink.targetID);
                if (surl.getPath().contains("/uniprot/") && turl.getPath().contains("/uniprot/")) {
                    String sid = surl.getPath().replace("/uniprot/", "");
                    String tid = turl.getPath().replace("/uniprot/", "");
                    String source = String.join(",", curLink.getDataSources());
                    String pubmed = String.join(",", curLink.getPublicationIDs(true));

                    String pathways = String.join(",", curLink.getPathwayUris());

                    sifPW.printf("%s\t%s\t%s\t%s\t%s\t%s\n",
                            sid, curLink.type, tid, source, pubmed, pathways);
                    if (curLink.type == SIFEnum.IN_COMPLEX_WITH) {
                        for (String mstr : curLink.getMediatorIDs()) {
                            URL murl = new URL(mstr);
                            String mid = murl.getPath().replace("/pc12/", "");
                            outComplex.write(mid + "\t" + sid + "\n");
                            outComplex.write(mid + "\t" + tid + "\n");
                        }
                    }
                }
            }
        }
        outSIF.close();
        outComplex.close();
    }
}
