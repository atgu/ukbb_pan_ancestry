---
id: background
title: Background
---

This FAQ is meant to accompany our public release of the results of a large genome-wide analysis for a wide range of outcomes in a diverse group of individuals from the UK. It is intended to provide context around and describe some of the limitations of our analyses to a lay audience. It does not go into comprehensive detail around every phenotype, but instead wishes to highlight the overarching goals of the project, how analyses were conducted, and describe some potential confounders that could affect our results. The practice of producing these public-facing FAQ documents was first adopted by the Social Science Genetics Association Consortium (https://www.thessgac.org/faqs) and has become a common practice among several genomics researchers. In this FAQ, we have followed a similar structure to that of others that have been previously released. For a more detailed description of the analyses and technical details, please refer to the wiki on our [github page](https://github.com/atgu/ukbb_pan_ancestry/wiki). We consider this a living document that we will update as the need arises. If you have questions or concerns about this FAQ, please reach out to us at [ukb.diverse.gwas@gmail.com](mailto:ukb.diverse.gwas@gmail.com).

## Who conducted this study?

:::note
[A team of researchers](/team) from the Analytic & Translational Genetics Unit (ATGU) at Massachusetts General Hospital and the Broad Institute of MIT and Harvard performed the analysis in this study.
:::note

The data used in these analyses are from the [UK Biobank](https://www.ukbiobank.ac.uk/), a large-scale open database with hundreds of thousands of individuals’ genotype data paired to electronic health records and survey measures. The UK Biobank recruited 500,000 people aged between 40-69 years in 2006-2010 from across the country to take part in this project. They have undergone measures, provided blood, urine and saliva samples for future analysis, detailed information about themselves, and agreed to have their health followed. The researchers at ATGU and the Broad Institute were not involved in the design of the UK Biobank resource or recruitment of participants, but have analyzed the breadth of this powerful resource.

Throughout this work, we regularly sought engagement and feedback from researchers and communities to help direct and contextualize this research and to discuss actions that will allow the substantial benefits of our analyses to outweigh the risks. Please see the section “[What has been done to reduce the potential harms of this research?](implications#what-has-been-done-to-reduce-the-potential-harms-of-this-research)” to learn more about the specific individuals and groups who have been a part of this effort and for information on what we have done to reduce risks to disadvantaged groups.


## What are the group's overarching goals?

:::note
Our fundamental goal was to build a resource that facilitates access to genetic association results (also known as summary statistics) for as many phenotypes in as many diverse populations as possible, particularly those that have traditionally been underrepresented in prior genetics work and excluded in most analyses of the widely-used UK Biobank resource.
:::note

We believe that easy access to high-quality data on diverse populations will accelerate research that can improve the health of the global population and can contribute to closing disparities that exist in the world. These association results can be used to better understand the biology underlying a broad range of traits. The additional specific benefits of including underrepresented populations in this study are described in detail in the section “[What are the potential benefits of this research?](implications#what-can-be-done-with-the-results-of-this-research-what-are-the-potential-benefits-of-this-research)” and throughout this document. Substantial health and social disparities exist between the groups studied in this research. While these disparities are largely a direct result of environmental factors, we hope that our research will lead to further work that mitigates disparities, even though we do not directly study those disparities as part of this project.


## Why was this study done?

:::note
This project is an effort to increase diversity in genetics research and to make use of data that is traditionally left out of analysis.
:::note

Historically in genetics research, for technical and social reasons [described below](study-design#why-have-certain-individuals-been-excluded-in-previous-research), most prior work has been done only in populations with European ancestries. Data from participants with predominantly non-European ancestries were usually excluded from previous studies including many studies of the UK Biobank. (See “[Why have certain individuals been excluded in previous research?](study-design#why-have-certain-individuals-been-excluded-in-previous-research)”)  Because of this, results from these genetic studies may not apply as well to other groups. This further implies that applying these findings to a clinical setting has the potential to increase, rather than decrease, health disparities.

Additionally, genetic studies require a large amount of expertise and computing power. By making these results publicly available, we remove this barrier to researchers and provide results for many phenotypes while using a consistent analysis pipeline. We hope this will accelerate the pace and reliability of genetic discoveries and encourage future studies to include data from more diverse participants.


## What is ancestry? Is it the same as race or ethnicity?

:::note
The “ancestry” of a group of people is related to the set of ancestors from whom they inherited their genetic variants. It does not have natural boundaries and it is <em>not</em> the same as race or ethnicity.
:::note

The distinctions between a person’s race and ancestry are important. “Ancestry” is a statistical construct based on the genetic variants that an individual inherited from their ancestors. “Race” and “ethnicity” are social constructs and group people based on _perceived_ shared physical, geographical, cultural, language, religion, or other social characteristics, often in an inherently unequal manner. As a result, a person’s race and ethnicity can depend on the time and place that an individual lives. Similarly, an individual’s self-identified race or ethnicity may at times differ from the corresponding genetic ancestry assigned by statistical algorithms. Treating ancestry, ethnicity, and race as equivalent concepts is incorrect. **In all our analyses, we exclusively refer to genetic ancestry.**

When measuring ancestry across individuals, geneticists use statistical tools and very large data to identify groups of people who are more genetically similar. Based on the region where people in a group live and what is known about human migration, populations with similar ancestries are often given geographic labels. For example, the vast majority of the individuals in the UK Biobank appear to have similar ancestries to other individuals who have grandparents native to countries in Northern Europe. For this reason, geneticists often refer to such individuals as having “European” ancestry. That said, **ancestry is a continuum that does not have obvious boundaries**. It is possible to divide a group of individuals into any number of “ancestry groups.” In this research, we use rough continental categories, as described in more detail in “_[How did you decide what populations to include? How did you assign individuals to each ancestry group?](study-design#how-did-you-decide-what-ancestry-groups-to-include-how-did-you-assign-individuals-to-each-ancestry-group)_”


## Terminology

### What is a GWAS?

A Genome-Wide Association Study (GWAS) is a scan of millions of genetic variants in the human genome, looking for variants that are associated with a particular “phenotype”.

### What is a phenotype?

A phenotype is a disease, outcome, or trait that can be measured and studied in a quantitative manner. Examples of phenotypes include height, self-reported smoking behavior, or whether a person has been diagnosed with type 2 diabetes.

### What is meant by "association"?

To test whether a genetic variant is associated with the phenotype, we compare individuals that have a copy or copies of the variant to those who have none. If the difference is large enough that it is very unlikely to have occurred by chance, we say that the variant is “associated” with the phenotype.


## What does it mean for a variant to be associated with a phenotype? Are the genetic variants discovered by GWAS “causal”?

:::note
Associations can occur for many reasons. Many of these reasons are not causal effects of the genetic variant on the phenotype.
:::note

For example, consider a case when the phenotype is a disease. It could be that the variant triggers a biological mechanism that directly leads to disease. This variant would be associated with the phenotype and would also be considered causal. On the other hand, it could be that the variant is by chance more common in certain groups or communities, and that they have a higher prevalence of disease due to environmental factors like pollution or social programs. This variant would also be associated with the phenotype but it would _not_ be considered causal. As a third example, many variants that have significant associations in a GWAS are not expected to themselves be causal, but instead might be associated because they are nearby on the chromosome to a variant that is causal. Finally, in some cases, genetic variants may only be associated with the phenotype in certain settings. Generally, GWAS can only tell us if a genetic variant is associated with a phenotype; **it cannot tell us why the variant is associated and does not demonstrate that the relationship is causal.**


## Do these results imply that genetics are responsible for the phenotypic differences between ancestry groups?

:::note
They do not. It is important to keep in mind that human populations are more genetically alike than we are different.
:::note

We don’t know whether there is a genetic basis for differences between groups because GWAS variants aren’t necessarily causal. Because of this, interpreting differences in genetic associations between groups is incredibly complicated, as associations can appear for a variety of reasons that don’t imply causality. Please see our explanation _[“What does it mean for a variant to be associated with a phenotype? Are the genetic variants discovered by GWAS “causal”?](background#what-does-it-mean-for-a-variant-to-be-associated-with-a-phenotype-are-the-genetic-variants-discovered-by-gwas-causal)_.

Additionally, for all of the phenotypes considered by this study,** there is not a single genetic variant that determines whether you will have a certain condition**. Instead, in most cases, there are many, many genetic variants that each have a small association with the phenotype, yielding on average very similar outcomes in different populations. Because these associations are so small, it takes substantial follow-up work and an enormous amount of data to determine whether a genetic variant has a population-specific effect -- often more data than are currently available to researchers.

This is further complicated by the fact that there are many other factors (unrelated to biological differences) that can lead to differences in associations between populations. Perhaps the most important of these reasons is that different populations often face different environmental factors (e.g., discrimination, rates of poverty, geography), and these factors may affect genetic associations. Together, this means that GWAS alone are insufficient to explain biological differences among populations.


## Since biology is mostly shared, why is diversity in genetics so important?

:::note
Diversity in genetic studies is important to ensure that findings are generalizable and that everyone can benefit from these findings.
:::note

While previous genetic studies have provided deep insights into the molecular basis of many phenotypes, participants in those studies have mostly consisted of groups of people who trace their ancestry back to regions within Europe. While all people have more genetic similarities than differences, certain genetic variants or combinations of variants are more common in groups that have ancestries that trace back to close-by regions. As a result, most previous studies have been best-suited for understanding the role of genetic variants that are more common in people from those European regions. Expanding genetics research to include individuals with diverse ancestries will improve our understanding of these phenotypes for everyone. For example, more diversity will help researchers establish which genetic variants are actually causal and which are just simply associated. It will also help researchers discover new biological mechanisms since some genetic variants are only common enough to study in certain populations. Discovering these biological mechanisms will help us better understand the underlying biology of important phenotypes shared by everyone. Additionally, studying underrepresented groups may make precision medicine possible for these currently underserved populations themselves.


## What did you learn as part of this study?

:::note
Our results from this study provide a starting point for understanding the genetic underpinnings of thousands of traits and diseases across global populations.
:::note

This study is ongoing and we are still actively analyzing results, but we intend to release these results well in advance of scientific publication or submission timelines because of the disproportionate value of these results to the field. For example, as part of our study, we are releasing GWAS results for a set of phenotypes related to COVID-19. While it is well-understood that differential rates of transmission of COVID-19 are primarily social and environmental, GWAS results may improve our understanding of the biological mechanisms that influence disease severity and accelerate the discovery of an urgently-needed treatment. By facilitating access to analyses conducted in underrepresented populations, we hope that the broader scientific community makes use of these data such that data from diverse ancestry participants in other cohorts are not discarded. As this study moves forward, we and others in the research community will continue to produce important scientific insights.
