---
id: study-design
title: Study Design
---

## What was done?

:::note
We used genetic and phenotypic data from ~500,000 participants in the UK Biobank to conduct [genome-wide association studies](background#what-is-a-gwas) (GWAS) among individuals with diverse ancestral backgrounds. This included more than 20,000 individuals with primarily non-European ancestries.
:::note

The [UK Biobank](https://www.ukbiobank.ac.uk/) is an open access database with hundreds of thousands of individuals’ genetic data paired to electronic health records and survey measures. We conducted GWAS for all phenotypes deemed to have sufficient statistical power. These phenotypes include a total of more than 16,000 GWAS conducted across a very broad range of phenotypes. A few examples of phenotypes we studied include anthropometric measures and physical attributes (e.g. height, BMI, bone density), blood panel traits (e.g. white blood cell count, cholesterol, blood glucose), common diseases (e.g. diabetes, cardiovascular disease, psychiatric disorders), electronic health record data (e.g. diagnosis codes entered by clinicians), prescription data (e.g. prescribed to take statins), health surveys (e.g. dietary intake, activity levels, general health satisfaction), social surveys (e.g. educational attainment, occupation), and many other measures. To summarize, phenotypes included both data pulled from electronic medical records as well as participants' survey responses to questionnaires given online or at the clinic.


## What data were used?

:::note
The data used in these analyses comes from the [UK Biobank](https://www.ukbiobank.ac.uk/), a large-scale open database with hundreds of thousands of individuals’ genotype data paired to electronic health records and survey data.
:::note

Researchers can gain access to the UK Biobank by writing a proposal for a research project, which then is reviewed for approval. The UK biobank is more thoroughly described on their website and in [a scientific publication](https://www.nature.com/articles/s41586-018-0579-z).

## Have you used data from countries other than the UK?

:::note
No, all participants lived in the UK at the time that the data were collected.
:::note

Because all of our data come from the UK Biobank and we did not collect additional data as part of this project, all individuals must have lived in the UK during the UK Biobank recruitment phase. (See “[How were participants recruited?](study-design#how-were-participants-recruited)”) However, the data do include individuals who were born outside of the UK but were recruited to be part of the UK Biobank sample after having moved to the UK. Since there are important differences in genetic diversity between the UK, the US, and other countries, we hope that resources comparable to this one will be produced as large new datasets become available in the future.


## How were participants recruited?

:::note
Our team did not recruit participants for this study. Instead, we analysed existing data from the [UK Biobank](https://www.ukbiobank.ac.uk/), a collection of ~500,000 individuals collected in the United Kingdom about ten years ago.
:::note

Here, we describe how participants in the UK Biobank cohort were recruited. Following the success of an initial pilot study in 2005-2006, the main stage of recruitment for the UK Biobank resource began in 2007, with the goal of recruiting 500,000 individuals between the ages of 40 and 69. The age restriction was due to the primary aim of the study: to improve the prevention, diagnosis and treatment of serious illnesses that typically onset later in life, including diabetes, cancer, arthritis, heart disease, stroke, and dementia. To that end, individuals from across the UK were contacted by post to participate in the study, with names, addresses, and dates of birth provided by the UK National Health Service (NHS). The 500,000 recruitment goal was reached in July of 2010, and recruitment ended shortly after. Focusing on voluntary recruitment of an older subset of the UK population sent by mail resulted in a sample of individuals that is more healthy and wealthy than the average Brit and that **had a greater fraction of European ancestries than the UK population overall**. This means that this cohort is not a perfect representation of the UK population, which further means that the results of our study may not reflect the UK population as a whole. This limitation is important to keep in mind when considering our results.


## How did you decide which phenotypes to study?

:::note
We included all phenotypes that were available to us for which there was sufficient data to conduct a GWAS.
:::note

Due to the scale of this project, we relied on quality control procedures that worked well in general for all phenotypes rather than using specific quality controls for each one. Quality controls are the procedures used to minimize errors in the data we use in our analysis. For example, individuals who report extremely large or small values of a phenotype are often just incorrectly recorded and can bias analyses. Therefore, we adjusted the values for such individuals using a standard transformation. For complete details on our quality control procedures, see our [wiki](https://github.com/atgu/ukbb_pan_ancestry/wiki). Certain phenotypes may require different quality control procedures than those we used in our analyses. As a result, some researchers may prefer to include different sets of individuals or define phenotypes slightly differently in their future work.


## How did you decide what ancestry groups to include? How did you assign individuals to each ancestry group?

:::note
The ancestry groups we used in these analyses are based on those described in two large existing globally and genetically diverse datasets. To assign each person to each ancestry group, we applied statistical methods to each individual’s genetic data.
:::note

Specifically, we compared the genome of each participant in the UK Biobank to the data in two large reference datasets containing genetically diverse individuals from across the globe, the [1000 Genomes Project](https://www.internationalgenome.org/) and [Human Genome Diversity Project (HGDP)](https://www.hagsc.org/hgdp/). We statistically assessed how similar each participant’s ancestry was to individuals from the populations included in these reference panels. These previous studies included labels which break down populations broadly into continental groupings that share ancestors and history over the course of tens to thousands of generations. The ancestry labels include African, American (which in these studies is the ancestry shared by many Hispanic/Latino groups), Central/South Asian, East Asian, European, and Middle Eastern. We assigned each individual into the ancestry groups that he/she was most similar to, adopting the same labels as used previously for consistency. We dropped those individuals who did not have a confident ancestry assignment. Notably, this approach does not rely on any other information, including self-reported race, ethnicity, or ancestry. We conducted our studies in all of the populations that had large enough numbers of individuals to learn about the genetic underpinnings of some traits, with individuals from each population analyzed together.

These ancestry labels have many limitations. First, as described in response to the question [“What is ancestry?”](background#what-is-ancestry-is-it-the-same-as-race-or-ethnicity), there is substantial diversity within each labeled population. Second, the 1000 Genome Project and HGDP data used to define continental populations have gaps in representation. They include more individuals from some regions than others, but this is not necessarily reflective of a region’s corresponding genetic diversity. For example, the individuals in the “African” population in the 1000 Genomes Project data have more West African ancestors than ancestors from other parts of Africa. For this reason, groups with ancestries from other parts of Africa may not be identified as accurately. Third, many individuals have ancestors from more than one of the groups defined by the 1000 Genome Project and HGDP. Such individuals are said to be “admixed,” which means that different parts of their genomes come from different continental populations. We discuss how such individuals affect our analysis more in the question below in [“What about people with mixed ancestries?”](study-design#what-about-people-with-mixed-ancestries). We hope that in the future, more diverse reference data sets will become available so that analyses such as ours will be less susceptible to these limitations.


## What about people with mixed ancestries?

:::note
Many admixed individuals are included in this study.
:::note

Since ancestry is a continuum (see _[“What is ancestry?”](background#what-is-ancestry-is-it-the-same-as-race-or-ethnicity)_), some ancestry labels consist nearly exclusively of individuals with varying amounts of admixture. For example, most individuals that are labeled as having “American” or “African” ancestries by our statistical algorithms share recent common ancestors with those labeled as having “European” ancestries. As long as there are enough people with a similar pattern of admixture -- as there are in the “African” and “American” ancestry groups -- we can study them together in a GWAS. However, some individuals have less common patterns of admixture, such that there are not enough similar individuals that we can group them together in a genetic study. Therefore, we currently drop such individuals from analysis. However, we believe that including these individuals should also be a priority, and plan to implement tools that allow us to include them in ongoing and future work to allow for increased inclusion of admixed participants in future studies.


## Why do you analyze ancestry groups separately?

:::note
Different ancestry groups are separately analyzed for statistical reasons; this does not imply that there are biological differences between ancestry groups (see [“Do these results imply that genetics are responsible for the phenotypic differences between ancestry groups?”](background#do-these-results-imply-that-genetics-are-responsible-for-the-phenotypic-differences-between-ancestry-groups)).
:::note

To help understand why previous scientific efforts have restricted to only one ancestry group, a classic GWAS example is of chopstick use. If we conducted such a GWAS in a group with some people with East Asian ancestries and some people with non-East Asian ancestries, we would almost surely find many associations. These associations would not likely correspond to biological mechanisms that affect manual dexterity or a personal preference for wooden cutlery, but they would just identify genetic variants that are _by chance_ more or less common in East Asian populations relative to the ancestries represented by others in the data. When GWAS is conducted in groups of individuals with similar ancestries, associations are less likely to be driven by these types of non-causal factors, which makes the results easier to interpret and makes follow-up work more productive.


## Why have certain individuals been excluded in previous research?

:::note
Due to the statistical reasons [described above](study-design#why-do-you-analyze-ancestry-groups-separately), there are some scientific advantages to studying large groups of genetically similar individuals. Unfortunately, due to a number of historical, infrastructural, political, ideological, monetary and other societal factors, this has resulted in the disproportionate recruitment of individuals with European ancestries, effectively excluding other groups from participating in most genetic studies.
:::note

[Previous work](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6494470/pdf/nihms-1506281.pdf) described some of the more societal explanations for why data collection has heavily focused on European populations. For statistical reasons (see ["Why do you analyze ancestry groups separately?"](study-design#why-do-you-analyze-ancestry-groups-separately)), standard practice in genetic studies has been to only analyze the largest homogenous subset of the data, which in practice means only including individuals who are of the largest ancestry group. Due to this Euro-centric bias in data collection, the largest group is most frequently those who are classified as having European ancestry. This snowballing effect has driven the overrepresentation of European-ancestry individuals in published GWAS, and perpetuates the continued data collection and study of European-ancestry individuals given that they have been more fully characterized.

A description of how ancestry is classified and where ancestry labels come from is in the section on [“How did you decide what populations to include? How did you assign individuals to each ancestry group?”](study-design#how-did-you-decide-what-ancestry-groups-to-include-how-did-you-assign-individuals-to-each-ancestry-group).


## Why are you including them now?

:::note
The UK Biobank provides a unique opportunity to study people with diverse ancestries. Although people with recent ancestors from outside the UK make up a small fraction of the data in the biobank, there are enough to run and learn new biology from GWAS. Many of these genotypic and phenotypic datasets are the largest available that include individuals with certain ancestries.
:::note

Because associations between single genetic variants and phenotypes tend to be very small, a large number of individuals need to be studied to find reliable associations. That is, we often need at least tens of thousands of individuals (and often more) to find and validate genetic associations. With this in mind, the groups of individuals that have been omitted from previous genetic studies were far too small to be able to produce reliable results. **Analyzing these data and releasing the results to the research community and the public will hopefully accelerate research that will benefit global populations.**

