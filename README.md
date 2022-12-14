# Evaluation of the robustness to batch effect in RNA velocity analyses

## Abstract

All analyses based on NGS data suffer from batch effect when there is more than one sample involved, for various reasons, such as different technologies used or carrying out the analysis at different times. We performed a series of simulations from synthetic and real datasets to study the effect of the batch effect on the analysis and conclusions from RNA velocity analyses with scVelo. The results of these simulations reveal that when the samples are correctly integrated and scVelo uses the principal components of the corrected counts in its spliced and unspliced count imputation algorithms, the conclusions of the analysis are not affected by the batch effect. When the batch effect is not confounded with any biological effect, RNA velocity analyses are suitable for integrated data, being able to take advantage of the extra power that comes of having more data.
