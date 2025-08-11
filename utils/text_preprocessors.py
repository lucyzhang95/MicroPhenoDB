import re


class TextStructurePreprocessor:
    """Rules for cleaning the structure of text strings."""

    def remove_special_char(self, name: str) -> str:
        """Removes special characters from the name."""
        return re.sub(r"[?!#*&+]", "", name).strip()

    def remove_colon4name(self, name: str) -> str:
        """Replaces colons with spaces in the name."""
        return re.sub(r":", " ", name).strip()

    def remove_dot4name_except_in_sp(self, name: str) -> str:
        """Removes dots from the name except in 'sp.' or 'spp.'."""
        name = re.sub(r"\b(sp|spp)\.", r"\1__dot__", name)
        numeric_matches = re.findall(r"\d+(?:\.\d+)+", name)
        for match in numeric_matches:
            protected = match.replace(".", "__dot__")
            name = name.replace(match, protected)
        name = name.replace(".", "")
        return name.replace("__dot__", ".").strip()

    def remove_hyphen4name(self, name: str) -> str:
        """Replaces hyphens with spaces in the name, handling specific cases."""
        name = name.replace("butyrate-producing", "BUTYRATEPRODUCING")
        name = re.sub(r"(?<=[a-z])-(?=\d)", " ", name)
        name = re.sub(r"(?<=[a-z])-(?=(like|associated|related|positive|negative|)\b)", " ", name)
        return name.replace("BUTYRATEPRODUCING", "butyrate-producing").strip()

    def split_name_by_slash(self, name: str) -> str:
        """Splits the name by slashes, keeping the first part."""
        return re.split(r"(?<=[a-zA-Z])/\s*(?=[a-zA-Z])", name)[0].strip()

    def split_on_conjunction_in_name(self, name, keyword, prefer):
        """Splits the name on a conjunction, returning the preferred part."""
        pattern = r"\b" + re.escape(keyword) + r"\b"
        parts = [p.strip() for p in re.split(pattern, name, maxsplit=1) if p.strip()]
        if not parts:
            return name.strip()
        if prefer == "left":
            return parts[0]
        if prefer == "right":
            return parts[1] if len(parts) > 1 else parts[0]
        return name

    def remove_parentheses(self, name: str) -> str:
        """Removes parentheses and their contents from the name."""
        return re.sub(r"\s*\(.+\)", "", name).strip()

    def process_comma(self, name: str) -> str:
        """Processes commas in the name, removing trailing words."""
        name = re.sub(r",\s*(and\s+)?[a-z]{1,3}\b", "", name)
        return name.split(",")[0].strip()


class TextSemanticPreprocessor:
    """Rules for cleaning the semantics of text strings."""

    def remove_non_english_chars_in_name(self, name: str) -> str:
        """Removes non-ASCII characters from the name."""
        return re.sub(r"[^\x00-\x7F]+", "", name)

    def remove_and_in_name(self, name: str) -> str:
        """Removes 'and' from the name, keeping the first part."""
        return name.split(" and ")[0].strip()

    def remove_in_and_one_word_after_in_name(self, name: str) -> str:
        """Removes 'in' followed by a single word from the name."""
        name = re.sub(r"\bin\b\s+\w+\b", "", name)
        return re.sub(r"\s{2,}", " ", name).strip()

    def remove_word_related_in_name(self, name: str) -> str:
        """Removes 'related' and similar terms from the name."""
        return re.sub(r"\b\w+-related\s*", "", name).strip()

    def remove_leading_non_in_name(self, name: str) -> str:
        """Removes leading non-taxonomic terms from the name."""
        return re.sub(r"^(?:non\w+\s+)+", "", name).strip()

    def remove_spp_in_name(self, name: str) -> str:
        """Removes 'sp.' or 'spp.' from the name."""
        return re.sub(r"\bspps?\b.*", "", name).strip()

    def remove_strain_in_taxon_name(self, name: str) -> str:
        """Removes 'strain' or 'strains' from the name."""
        return re.sub(r"\s{2,}", " ", re.sub(r"\bstrains?\b", "", name).strip())

    def remove_type_in_taxon_name(self, name: str) -> str:
        """Removes 'type' or 'types' from the name."""
        return re.sub(r"\btypes?\b", "", name).strip()

    def remove_group_in_taxon_name(self, name: str) -> str:
        """Removes 'group' or 'groups' from the name."""
        name = re.sub(r"\bgroups?\s+[a-z]\b", "", name)
        name = re.sub(r"\b(groups?|subgroup)\b$", "", name)
        return re.sub(r"\b(groups?|subgroup)\b(?=\s+[a-z])", "", name).strip()

    def remove_serovar_in_taxon_name(self, name: str) -> str:
        """Removes 'serovar' or 'serovars' from the name."""
        return re.sub(r"\bserovars?.*\b", "", name).strip()

    def remove_pre_postfix_in_taxon_name(self, name: str) -> str:
        """Removes common prefixes and suffixes from the name."""
        name = re.sub(r"^\s*b\s+", "", name)
        name = re.sub(r"\bstains?\b", "", name)
        name = re.sub(
            r"(coagulase negative|(?:non)?hemolytic|sensu lato|complex(?:es)?|incertae sedis|rapid growers)",
            "",
            name,
        )
        return re.sub(r"\s{2,}", " ", name).strip()

    def expand_taxon_name_abbrev(self, name: str) -> str:
        """Expands common abbreviations in taxon names."""
        expansions = {
            r"\be histolytica\b": "entamoeba histolytica",
            r"\bp multocida\b": "pasteurella multocida",
            r"\bhsv\b": "herpes simplex virus",
            r"\bhpv\b": "human papillomavirus",
            r"\bhiv\b": "human immunodeficiency virus",
            r"\blcm\b": "lymphocytic choriomeningitis",
            r"\bcluster\b": "family",
            r"\bspirochaeta\b": "borrelia",
            r"\bpiv\b": "parainfluenza virus",
            r"\btm7\b": "candidatus saccharimonadota",
            r"\brubella\b": "rubella virus",
            r"\bmumps\b": "mumps virus",
            r"\bntm\b": "mycobacteriales",  # changed from "non-tuberculous mycobacteria" to "mycobacteriales"
            r"^\bsr1\b$": "candidatus absconditibacteriota",
            r"zygomycetes": "mucoromycota",  # zygomycetes is an obsolete term for mucoromycota and zoopagomycota.
        }
        for pattern, replacement in expansions.items():
            name = re.sub(pattern, replacement, name)
        return re.sub(r"\s{2,}", " ", name).strip()

    def replace_taxon_name(self, name: str) -> str:
        """Standardizes plural or variant taxon names to their singular form and correct spelling."""
        replacements = {
            r"streptococci": "streptococcus",
            r"lactobacilli": "lactobacillus",
            r"enterococci": "enterococcus",
            r"staphylococci": "staphylococcus",
            r"coxsackie": "coxsackievirus",
            r"gemellales": "gemella",
            r"\bparainfluenza virus\b": "orthorubulavirus",
            r"\bparainfluenza\b": "orthorubulavirus",  # make rank broader to include all parainfluenza viruses
            r"\bpapovavirus\b": "papillomavirus",  # papovavirus obsolete term for both papillomavirus and polyomavirus
        }
        for old, new in replacements.items():
            name = re.sub(rf"\b{old}\b", new, name, flags=re.IGNORECASE)
        return name.strip()
