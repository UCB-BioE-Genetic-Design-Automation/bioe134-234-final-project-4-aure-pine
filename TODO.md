- Adding in Cas9 restriction enzyme functionality?
Adding other cloning techniques?

### Functions:
- `find_restriction_sites`
- `combine_sequences`
    - to combine input CDS, UTRs etc into a full sequence
- `annotate_seq`
    - Using biopython methods: `SeqFeature` and `SeqRecord`
        - Create features for:
            - Promoters
            - UTRs (5’ or 3’)
            - Terminators
            - Restriction sites


1. Check the component names, if its not one of the list, print a warning that they will be given the feature name: misc_feature
2. Combine the sequences together: CDS, with the components