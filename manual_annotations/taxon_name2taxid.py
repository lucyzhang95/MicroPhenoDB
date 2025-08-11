"""
This file contains a mapping of taxon names to their corresponding NCBI Taxonomy IDs (taxids).
"""

MANUAL_MAP_UNMAPPED_TAXON_NAMES = {
    "actinomyces neuii subspecies neuii": 144053,
    "actinomyces radingae": 131110,
    "actinomyces turicensis": 131111,
    "candidate division candidatus saccharimonadota genomosp": 239137,
    "candidate division candidatus saccharimonadota single cell isolate tm7b": 447455,
    "candidate division candidatus saccharimonadota single cell isolate tm7c": 447456,
    "clostridia family i": 31979,
    "clostridium family xiva": 543317,
    "clostridium family xviii": 189325,
    "clostridium lituseburense": 1537,
    "clostridium xi": 186804,
    "clostridium xivb": 543317,
    "coxsackievirus a virus": 12066,
    "creutzfeldt jakob disease": 36469,
    "cryptococcus albidus": 100951,
    "cysticercosis": 6204,
    "escherichia vulneris": 566,
    "eubacterium tortuosum": 39494,
    "hookworms cestodes": 6157,
    "neisseriagonorrhoeae": 485,
    "orthorubulavirus 1c4": 2560526,
    "ovine jaagziekte virus": 11746,
    "prevotella multisaccharivorax": 310514,
    "prevotella nanceiensis": 425941,
    "saccharomyces castellii": 27288,
    "syphilis": 160,
    "trophyrema": 2039,
    "uncultured clostridiales ii": 186801,
    "vulvovaginal candidiasis": 5476,
}

OVERRIDE_BT_mapped_TAXID = {
    "influenza": {"taxid": 11320},  # influenza a virus
    "bifidobacterium infantis": {"taxid": 1682},
    "powassan": {"taxid": 11083},  # powassan virus
    "rubella virus virus": {"taxid": 11041},  # rubella virus
    "st louis encephalitis": {"taxid": 11080},  # st. louis encephalitis virus
    "yeasts": {"taxid": 5206},  # cryptococcus
}

