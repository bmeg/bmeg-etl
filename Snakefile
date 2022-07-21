


rule all:
	input:
		"source/pathway_commons/PathwayCommons12.bind.complex",
		"source/pathway_commons/PathwayCommons12.humancyc.extSIF",
		"source/pathway_commons/PathwayCommons12.humancyc.complex",
		"source/pathway_commons/PathwayCommons12.dip.extSIF",
		"source/pathway_commons/PathwayCommons12.pid.complex",
		"source/pathway_commons/PathwayCommons12.msigdb.extSIF",
		"source/pathway_commons/PathwayCommons12.kegg.complex",
		"source/pathway_commons/PathwayCommons12.innatedb.complex",
		"source/pathway_commons/PathwayCommons12.pid.extSIF",
		"source/pathway_commons/PathwayCommons12.msigdb.complex",
		"source/pathway_commons/PathwayCommons12.innatedb.extSIF",
		"source/pathway_commons/PathwayCommons12.dip.complex",
		"source/pathway_commons/PathwayCommons12.inoh.complex",
		"source/pathway_commons/PathwayCommons12.corum.extSIF",
		"source/pathway_commons/PathwayCommons12.psp.complex",
		"source/pathway_commons/PathwayCommons12.mirtarbase.extSIF",
		"source/pathway_commons/PathwayCommons12.hprd.complex",
		"source/pathway_commons/PathwayCommons12.mirtarbase.complex",
		"source/pathway_commons/PathwayCommons12.pathbank.extSIF",
		"source/pathway_commons/PathwayCommons12.pathbank.complex",
		"source/pathway_commons/PathwayCommons12.psp.extSIF",
		"source/pathway_commons/PathwayCommons12.biogrid.complex",
		"source/pathway_commons/PathwayCommons12.reconx.complex",
		"source/pathway_commons/PathwayCommons12.corum.complex",
		"source/pathway_commons/PathwayCommons12.drugbank.complex",
		"source/pathway_commons/PathwayCommons12.bind.extSIF",
		"source/pathway_commons/PathwayCommons12.reconx.extSIF",
		"source/pathway_commons/PathwayCommons12.biogrid.extSIF",
		"source/pathway_commons/PathwayCommons12.panther.complex",
		"source/pathway_commons/PathwayCommons12.kegg.extSIF",
		"source/pathway_commons/PathwayCommons12.reactome.complex",
		"source/pathway_commons/PathwayCommons12.ctd.extSIF",
		"source/pathway_commons/PathwayCommons12.panther.extSIF",
		"source/pathway_commons/PathwayCommons12.drugbank.extSIF",
		"source/pathway_commons/PathwayCommons12.inoh.extSIF",
		"source/pathway_commons/PathwayCommons12.netpath.extSIF",
		"source/pathway_commons/PathwayCommons12.reactome.extSIF",
		"source/pathway_commons/PathwayCommons12.netpath.complex",
		"source/pathway_commons/PathwayCommons12.ctd.complex",
		"source/pathway_commons/PathwayCommons12.hprd.extSIF"

rule :
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.msigdb.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.msigdb.extSIF",
		"source/pathway_commons/PathwayCommons12.msigdb.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.msigdb.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.msigdb"

rule _1:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.panther.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.panther.extSIF",
		"source/pathway_commons/PathwayCommons12.panther.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.panther.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.panther"

rule _2:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.pathbank.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.pathbank.extSIF",
		"source/pathway_commons/PathwayCommons12.pathbank.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.pathbank.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.pathbank"

rule _3:
	output:
		"source/pathway_commons/PathwayCommons12.dip.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.dip.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.dip.BIOPAX.owl.gz"

rule _4:
	output:
		"source/pathway_commons/PathwayCommons12.humancyc.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.humancyc.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.humancyc.BIOPAX.owl.gz"

rule _5:
	output:
		"source/pathway_commons/PathwayCommons12.panther.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.panther.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.panther.BIOPAX.owl.gz"

rule _6:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.kegg.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.kegg.extSIF",
		"source/pathway_commons/PathwayCommons12.kegg.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.kegg.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.kegg"

rule _7:
	output:
		"source/pathway_commons/PathwayCommons12.innatedb.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.innatedb.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.innatedb.BIOPAX.owl.gz"

rule _8:
	output:
		"source/pathway_commons/PathwayCommons12.bind.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.bind.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.bind.BIOPAX.owl.gz"

rule _9:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.bind.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.bind.extSIF",
		"source/pathway_commons/PathwayCommons12.bind.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.bind.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.bind"

rule _10:
	output:
		"source/pathway_commons/PathwayCommons12.pathbank.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.pathbank.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.pathbank.BIOPAX.owl.gz"

rule _11:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.innatedb.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.innatedb.extSIF",
		"source/pathway_commons/PathwayCommons12.innatedb.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.innatedb.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.innatedb"

rule _12:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.dip.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.dip.extSIF",
		"source/pathway_commons/PathwayCommons12.dip.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.dip.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.dip"

rule _13:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.reconx.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.reconx.extSIF",
		"source/pathway_commons/PathwayCommons12.reconx.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.reconx.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.reconx"

rule _14:
	output:
		"source/pathway_commons/PathwayCommons12.hprd.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.hprd.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.hprd.BIOPAX.owl.gz"

rule _15:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.corum.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.corum.extSIF",
		"source/pathway_commons/PathwayCommons12.corum.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.corum.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.corum"

rule _16:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.hprd.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.hprd.extSIF",
		"source/pathway_commons/PathwayCommons12.hprd.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.hprd.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.hprd"

rule _17:
	output:
		"source/pathway_commons/PathwayCommons12.corum.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.corum.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.corum.BIOPAX.owl.gz"

rule _18:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.drugbank.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.drugbank.extSIF",
		"source/pathway_commons/PathwayCommons12.drugbank.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.drugbank.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.drugbank"

rule _19:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.reactome.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.reactome.extSIF",
		"source/pathway_commons/PathwayCommons12.reactome.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.reactome.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.reactome"

rule _20:
	output:
		"source/pathway_commons/PathwayCommons12.msigdb.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.msigdb.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.msigdb.BIOPAX.owl.gz"

rule _21:
	output:
		"source/pathway_commons/PathwayCommons12.biogrid.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.biogrid.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.biogrid.BIOPAX.owl.gz"

rule _22:
	output:
		"source/pathway_commons/PathwayCommons12.kegg.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.kegg.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.kegg.BIOPAX.owl.gz"

rule _23:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.pid.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.pid.extSIF",
		"source/pathway_commons/PathwayCommons12.pid.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.pid.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.pid"

rule _24:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.psp.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.psp.extSIF",
		"source/pathway_commons/PathwayCommons12.psp.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.psp.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.psp"

rule _25:
	output:
		"source/pathway_commons/PathwayCommons12.inoh.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.inoh.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.inoh.BIOPAX.owl.gz"

rule _26:
	output:
		"source/pathway_commons/PathwayCommons12.netpath.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.netpath.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.netpath.BIOPAX.owl.gz"

rule _27:
	output:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar"
	shell:
		"cd transform/pathway_commons && docker run  -u `id -u` -v `pwd`:`pwd` -v `pwd`/maven-repo:/root/.m2 -w `pwd` -ti maven:3-openjdk-11 mvn package"

rule _28:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.mirtarbase.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.mirtarbase.extSIF",
		"source/pathway_commons/PathwayCommons12.mirtarbase.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.mirtarbase.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.mirtarbase"

rule _29:
	output:
		"source/pathway_commons/PathwayCommons12.pid.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.pid.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.pid.BIOPAX.owl.gz"

rule _30:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.humancyc.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.humancyc.extSIF",
		"source/pathway_commons/PathwayCommons12.humancyc.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.humancyc.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.humancyc"

rule _31:
	output:
		"source/pathway_commons/PathwayCommons12.reconx.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.reconx.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.reconx.BIOPAX.owl.gz"

rule _32:
	output:
		"source/pathway_commons/PathwayCommons12.ctd.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.ctd.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.ctd.BIOPAX.owl.gz"

rule _33:
	output:
		"source/pathway_commons/PathwayCommons12.mirtarbase.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.mirtarbase.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.mirtarbase.BIOPAX.owl.gz"

rule _34:
	output:
		"source/pathway_commons/PathwayCommons12.reactome.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.reactome.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.reactome.BIOPAX.owl.gz"

rule _35:
	output:
		"source/pathway_commons/PathwayCommons12.drugbank.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.drugbank.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.drugbank.BIOPAX.owl.gz"

rule _36:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.inoh.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.inoh.extSIF",
		"source/pathway_commons/PathwayCommons12.inoh.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.inoh.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.inoh"

rule _37:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.netpath.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.netpath.extSIF",
		"source/pathway_commons/PathwayCommons12.netpath.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.netpath.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.netpath"

rule _38:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.ctd.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.ctd.extSIF",
		"source/pathway_commons/PathwayCommons12.ctd.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.ctd.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.ctd"

rule _39:
	input:
		"transform/pathway_commons/target/pc-extract-1.0-jar-with-dependencies.jar",
		"source/pathway_commons/PathwayCommons12.biogrid.BIOPAX.owl.gz"
	output:
		"source/pathway_commons/PathwayCommons12.biogrid.extSIF",
		"source/pathway_commons/PathwayCommons12.biogrid.complex"
	resources:
		mem_mb=15000
	shell:
		"cd transform/pathway_commons && java -Xmx15g -jar target/pc-extract-1.0-jar-with-dependencies.jar ../../source/pathway_commons/PathwayCommons12.biogrid.BIOPAX.owl.gz ../../source/pathway_commons/PathwayCommons12.biogrid"

rule _40:
	output:
		"source/pathway_commons/PathwayCommons12.psp.BIOPAX.owl.gz"
	shell:
		"cd transform/pathway_commons && curl -o ../../source/pathway_commons/PathwayCommons12.psp.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.psp.BIOPAX.owl.gz"



