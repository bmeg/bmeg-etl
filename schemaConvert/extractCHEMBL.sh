gunzip -c pharmacodb/pharmacodb_profiles_gx.drugResponse.drug_response.json.gz | jq ' .compounds | .[]?|.id' | sort | uniq | less -S
comm -12 <(gunzip -c g2p/g2p.main.assocation.json.gz | jq ' .compounds | .[]?|.id' | sort | uniq) <(gunzip -c pharmacodb/pharmacodb_profiles_gx.drugResponse.drug_response.json.gz | jq ' .compounds | .[]?|.id' | sort | uniq) | wc -l
