RegionJobTemplate = """cd %(location)s ; \
nice -n14 ~/rosetta/rosetta_features/rosetta.gcc \
%(series)s %(protein)s %(chain)s -paths %(paths)s \
-score \
-scorefile data/%(group)s_region \
-decoyfeatures region \
-fa_input \
-try_both_his_tautomers \
-l %(location)s/%(group)s.list 
"""

FeaturesJobTemplate = """cd %(location)s ; \
rm -f data/%(group)s_features.fasc ; \
nice -n14 ~/rosetta/rosetta_features/rosetta.gcc \
%(series)s %(protein)s %(chain)s -paths %(paths)s \
-score \
-scorefile data/%(group)s_features \
-decoyfeatures features \
-fa_input \
-try_both_his_tautomers \
-l %(location)s/%(group)s.list 
"""
