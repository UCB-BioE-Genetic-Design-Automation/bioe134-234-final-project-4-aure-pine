from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis, EcoRI, BamHI

seq = Seq("ATGCGCAGCTGTTCAGAGTTCATGAATTCTTATAATGGTAATTTCAAAGCAAGATTGGGAGCACAATTGTCAAATCTATGATGAATTCTATTTTGTTGCGTTATACCCAGCTTTTGGGCTTTTTCGAGAACGATCTCCCT")

rb = RestrictionBatch([EcoRI, BamHI])
anal = Analysis(rb, seq)
anal.print_as('map')
anal.print_that()

