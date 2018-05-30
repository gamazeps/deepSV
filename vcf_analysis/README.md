# VCF analysis for generating negative classes

The VCF file provides us with variants for our model, however if we want to be able to discriminate
between variants and non variants we need examples of regions that do not contain variants.

The choice of these regions has its importance as the model will consider them for learning
discrimination rules, we thus need to provide non variants that look close to what our variant
caller will encounter in real use.

I am currently not aware of an easy way to get the candidate variants out of variant callers, if we
were to be able to get them they would provide the best negative samples in my opinion.

## Characteristics of good negative classes

We would like our model not to just learn the region where variants happen, indeed the number of
regions is not that high, and our model could just end up learning them.

This would behave the same as a decision function that just looks if a candidate region was ever
included in a VCF file.

*We will thus need to generate tensors of regions where variants happen, in humans where that
variant is not present*, thus "un-biasing" the model for locus specificity.

We would also like to have regions where no variants have ever been found, as they would be closest
to "non variant" as imaginable.

We would also like to have regions from the previous two categories, with artifacts and noise
artificially introduced in the reads, in order to make our model more robust to them.

Ideally we would also like to get the false positives of other variant callers and train our model
to use them, this last part should be evaluated separately in order to avoid an unfair advantage for
our model in evaluation.
