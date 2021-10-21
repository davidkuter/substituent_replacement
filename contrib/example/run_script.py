from substituent_replacement.api.substituent_replacement import substituent_replacement


substituent_replacement(target_smiles="c1cc(ccc1c2c(n(cn2)CC3CC3)c4ccnc(n4)N)F",  # Pat Walters example
                        core_smarts="n1cncc1",                                    # Pat Walters example
                        out_folder='/home/dkuter/scratch/results',
                        level=1)                                                  # Options are 1 and 2
