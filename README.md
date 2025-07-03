# Data parser for [MicroPhenoDB](http://liwzlab.ifr.fidt.top:61010/microphenodb/#/home)
A manually collected and curated database includes **Microbe-Disease** associations (5529 from downloads) from HMDAD and Disbiome databases. 

## Statistics Summary
Last Modified: 05/14/2025 <br>
**Number of records: 5342, 187 records fewer than the full data (5529). <br>
~96.6% of the entire records** <br>
Some of the taxon and disease preprocessing workflow can refer to

<details>
<summary>Here</summary>

https://docs.google.com/spreadsheets/d/1dnPfB6qppecZWK3Yl_6HHXM1M4i7Vpam55-_vkHo2DI/edit?gid=0#gid=0

</details>

### Taxon 
- Total unique taxon names: 1767
- NCIT cached: 582 
- ETE3 cached: 1031 
- ENTREZ cached: 56 
- BT cached: 22 
- RapidFuzz cached: 53
- Unique taxids for query to get full taxon info: 1621
- Complete taxon mapped: 1813
- Taxon with description count: 2985


```json
{
   "genus":2345,
   "species":2315,
   "family":335,
   "no rank":118,
   "phylum":106,
   "serotype":40,
   "order":34,
   "class":18,
   "subspecies":11,
   "serogroup":9,
   "kingdom":4,
   "biotype":3,
   "strain":2,
   "forma":1,
   "subgenus":1
}
```

### Disease
- Total disease names: 5518
- All unique disease mapped: 485
- Disease with descriptions: 5336

```json
{
   "EFO":4680,
   "MONDO":301,
   "orphanet":249,
   "HP":108,
   "DOID":4
}
```

### Publication
- Unique pmids: 515
- Fetched pubmed metadata: 514

### Anatomical Entity
```json
{
   "gastrointestinal tract":2096,
   "other":1350,
   "oral cavity":541,
   "respiratory tract":250,
   "skin and soft tissue":222,
   "urinary tract":173,
   "vagina":141,
   "central nervous system":112,
   "bloodstream":80,
   "nasal cavity":75,
   "eye":63,
   "throat":56,
   "genitals":46,
   "cervix":30,
   "placenta":24,
   "gallbladder":22,
   "oesophagus":12,
   "foot":11,
   "unknown":11,
   "abdominal cavity":10,
   "auditory meatus":7,
   "biofilm":6,
   "lung":3,
   "breast":1
}
```


## Output Example
```json
{
   "_id":"496a78cb-c71f-4e78-a997-80ed0be4e241",
   "association":{
      "predicate":"biolink:OrganismalEntityAsAModelOfDiseaseAssociation",
      "type": "biolink:associated_with",
      "qualifier":"increase",
      "score":0.599701933,
      "anatomical_entity":"vagina",
      "infores":"MicroPhenoDB",
      "publication":{
         "id":"PMID:21272848",
         "pmid":21272848,
         "name":"Microarray-based identification of clinically relevant vaginal bacteria in relation to bacterial vaginosis.",
         "summary":"The objective was to examine the use of a tailor-made DNA microarray containing probes representing the vaginal microbiota to examine bacterial vaginosis. One hundred one women attending a health center for HIV testing in South Africa were enrolled. Stained, liquid-based cytology slides were scored for bacterial vaginosis. An inventory of organisms was obtained using microarray technology, probing genera associated with bacterial vaginosis in more detail, namely Gardnerella, Atopobium, Dialister, Leptotrichia, Megasphaera, Mobiluncus, Peptostreptococcus, Prevotella, and Sneathia. Of 101 women, 34 were diagnosed positive for bacterial vaginosis. This condition was associated with an increased microbial diversity. It is no longer useful to base the diagnosis of bacterial vaginosis on Gardnerella alone. Rather, its presence with Leptotrichia and Prevotella species, and especially Atopobium was more indicative of an aberrant state of the vaginal flora. To understand the vaginal microbiota in more detail, microarray-based identification can be used after microscopic scoring. [abstract]",
         "doi":"10.1016/j.ajog.2010.11.012",
         "type":"biolink:Publication"
      }
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
   }
}
```

## Issues and Future Updates
Most of the clostridium clusters are manually mapped to `taxid:189325`
*Clostridia incertae sedis*. For further mapping in depth with cluster or family along with the species or strains under the cluster requires more knowledge regarding the evolution and taxonomy of 
*Clostridia*. Can refer to [this paper]("https://pmc.ncbi.nlm.nih.gov/articles/PMC6656338/").

## Biothings Data plugin Type & Stats Report
The data plugin indicated a few duplicated _ids when I was using subject and object ids.
I have randomly checked some ids and they are from different publications, so I had to use uuid to generate random ids to prevent merging or replacement of the records.