# Data parser for [MicroPhenoDB](http://liwzlab.ifr.fidt.top:61010/microphenodb/#/home)
A manually collected and curated database includes microbe-disease (5529 from downloads) associations from HMDAD and Disbiome databases. 

## Example of output record
```ruby
{
   "_id":"2051_OrganismalEntityAsAModelOfDiseaseAssociation_0005316",
   "association":{
      "predicate":"biolink:OrganismalEntityAsAModelOfDiseaseAssociation",
      "qualifier":"increase",
      "score":0.599701933,
      "anatomical_entity":"vagina",
      "infores":"MicroPhenoDB"
   },
   "object":{
      "id":"MONDO:0005316",
      "name":"bacterial vaginosis",
      "description":"Infection caused by bacterial overgrowth in the vagina. Most affected women are asymptomatic. When symptoms occur, they include foul-smelling vaginal discharge, vaginal itching, and burning. Risk factors include sexual activity with multiple partners and the use of vaginal douches and intrauterine devices. Up to a third of cases resolve without treatment. Antibiotic treatment is recommended when symptoms are present and for women that are pregnant at the time of infection. [NCIT:P378]",
      "type":"biolink:Disease",
      "xrefs":{
         "mondo":"MONDO:0005316"
      },
      "original_name":"bacteria vaginosis"
   },
   "subject":{
      "id":"taxid:2051",
      "taxid":2051,
      "description":"A species of anaerobic, Gram positive, curved rod shaped bacterium assigned to the phylum Actinobacteria. This species is motile by one to six flagella that originate from the same spot on each cell and is oxidase, indole and catalase negative. M. curtisii is found in the vaginal tract and is pathogenic, being a causative agent of bacterial vaginosis.[NCIT]",
      "name":"mobiluncus curtisii",
      "parent_taxid":2050,
      "lineage":[
         2051,
         2050,
         2049,
         2037,
         1760,
         201174,
         1783272,
         2,
         131567,
         1
      ],
      "rank":"species",
      "xrefs":{
         "ncit":"C86518"
      },
      "type":"biolink:OrganismTaxon",
      "original_name":"mobiluncus curtisii"
   },
   "publication":{
      "pmid":21272848,
      "title":"Microarray-based identification of clinically relevant vaginal bacteria in relation to bacterial vaginosis.",
      "abstract":"The objective was to examine the use of a tailor-made DNA microarray containing probes representing the vaginal microbiota to examine bacterial vaginosis. One hundred one women attending a health center for HIV testing in South Africa were enrolled. Stained, liquid-based cytology slides were scored for bacterial vaginosis. An inventory of organisms was obtained using microarray technology, probing genera associated with bacterial vaginosis in more detail, namely Gardnerella, Atopobium, Dialister, Leptotrichia, Megasphaera, Mobiluncus, Peptostreptococcus, Prevotella, and Sneathia. Of 101 women, 34 were diagnosed positive for bacterial vaginosis. This condition was associated with an increased microbial diversity. It is no longer useful to base the diagnosis of bacterial vaginosis on Gardnerella alone. Rather, its presence with Leptotrichia and Prevotella species, and especially Atopobium was more indicative of an aberrant state of the vaginal flora. To understand the vaginal microbiota in more detail, microarray-based identification can be used after microscopic scoring.",
      "doi":"10.1016/j.ajog.2010.11.012",
      "type":"biolink:Publication"
   }
}
```
