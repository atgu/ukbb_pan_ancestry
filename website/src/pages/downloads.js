import React from 'react'
import classnames from 'classnames'
import Layout from '@theme/Layout'
import Link from '@docusaurus/Link'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import useBaseUrl from '@docusaurus/useBaseUrl'
import styles from './styles.module.css'

function Downloads() {
  const context = useDocusaurusContext()
  const { siteConfig = {} } = context
  return (
    <Layout title={`${siteConfig.title}`} description="Downloads">
      <header>
        <div className="container">
          <h1 className="page-title">Downloads</h1>
        </div>
      </header>
      <main>
        <div className="container">
          <p className="">
            The GWAS results are available for download in two main formats:
          </p>
          <ul>
            <li>
              <b>Per-phenotype flat files</b>: for most analyses of <b>one or a few phenotypes</b>, we suggest using the
              per-phenotype flat files, available freely on Amazon AWS. More information on the file formats is
              available in the <Link to={useBaseUrl('docs/per-phenotype-files')}>Technical Details</Link>.
              <ul>
                <li>
                  The phenotype manifest (browse on <a href="https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=903887429" target="_blank">Google Sheets</a> or
                  download on <a href="https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz">Amazon AWS</a>) contains the location and
                  detailed information of all per-phenotype files for those phenotypes for which GWAS was run.
                </li>
                <li>
                  The variant manifest contains detailed information on each variant
                  (download on <a href="https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz">Amazon AWS</a>
                  , <a href="https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz.tbi">tbi</a>).
                </li>
              </ul>
            </li>
            <li>
              <b>Hail format</b>: For large-scale analyses of many phenotypes, we provide the full dataset
              in <Link to={useBaseUrl('docs/hail-format')}>Hail MatrixTable format</Link> on Google Cloud.
            </li>
          </ul>
          <p className="">
            In addition, the LD matrices and scores are available in the following formats:
          </p>
          <ul>
            <li>
              <b>LDSC-compatible flat files</b>: for running LD score regression, we suggest using the
              LD score flat files available on Amazon AWS (download the tarball file <a href="https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz">here</a>). More information on the file formats is
              available on <a href="https://github.com/bulik/ldsc/wiki" target="_blank">the LDSC website</a>.
            </li>
            <li>
              <b>Hail format</b>: For large-scale analyses, we provide the full LD matrices and scores
              in <Link to={useBaseUrl('docs/hail-format')}>Hail format</Link> on Amazon AWS.
            </li>
          </ul>
          <p className="">
            All heritability estimates (see <Link to={useBaseUrl('docs/heritability')}>here</Link> for more information on our approach) are available for download in the following formats:
          </p>
          <ul>
            <li>
              <b>Flat files</b>: the manifest flat file is available on AWS (<a href="https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/h2_manifest.tsv.bgz">tarball here</a>) or on <a href="https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1797288938" target="_blank">Google Sheets</a>.
              Our topline results can be found as part of the main phenotype manifest (<a href="https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz">Amazon AWS</a> or <a href="https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=903887429" target="_blank">Google Sheets</a>)
            </li>
            <li>
              <b>Hail format</b>: for large-scale analyses and integration with our other datasets we provide heritability data
              in <Link to={useBaseUrl('docs/hail-format')}>Hail format</Link> on Google Cloud Platform.
            </li>
          </ul>
          <p className="">
            The phenotype correlation matrix, used in the construction of the maximally independent set of phenotypes passing QC (see <Link to={useBaseUrl('blog/2022/04/11/h2-qc-updated-sumstats')}>here</Link> for details),
            is available on Amazon AWS (download the tarball file <a href="https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/Pheno_pairwise_correlations_update.txt.bgz">here</a>).
          </p>
          <p className="">
            The ancestry assignments (as well as corresponding principal components and covariates used in our analyses)
            are available for download through the UK Biobank portal
            as <a href="https://biobank.ctsu.ox.ac.uk/showcase/dset.cgi?id=2442">Return 2442</a>.
            These are available to researchers registered with the UK Biobank:
            refer to instructions within the AMS portal to download these results.
          </p>
          <h3>Terms</h3>
          <p>
            All data here are released openly and publicly for the benefit of the wider biomedical community.
            You can freely download and search the data, and we encourage the use and publication of results
            generated from these data. From the perspective of the Pan-UKB team, there are absolutely no restrictions
            or embargoes on the publication of results derived from this data. However, we note that this research has
            been conducted using the UK Biobank Resource (project ID 31063), and use of this data is bound by all terms
            of usage of the UK Biobank: more information about can be
            found <a href="https://www.ukbiobank.ac.uk/principles-of-access">here</a>.
            All users of this data agree to not attempt to reidentify participants.
          </p>
          <p>
            These data are provided on an "AS-IS" basis, without warranty of any type, expressed or implied, including but not limited to any warranty as to their performance, merchantability, or fitness for any particular purpose (see license information below). This dataset has been subjected to quality control, but variant calling and statistical methods to associate variants and phenotypes is an imperfect and probabilistic process, so many errors no doubt remain: if you find any glaring errors, feel free to contact us. Users of the dataset certify that they are in compliance with all applicable local, state, and federal laws or regulations and institutional policies regarding human subjects and genetics research.
          </p>
          <p>
            The GWAS results data produced by the Pan-UKB are available free of restrictions under
            the <a href="https://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International
            (CC BY 4.0)</a>. We request that you acknowledge and give attribution to both the Pan-UKB project and
            UK Biobank, and link back to the relevant page, wherever possible.
          </p>
          <h3>Citation</h3>
          <p>
            In addition to acknowledging the UK Biobank, we request that any use of this dataset in publications cite:
          </p>
          <p><code>
            Pan-UKB team. https://pan.ukbb.broadinstitute.org. 2020.
          </code></p>
          <p>
            There is no need to include us as authors on your manuscript, unless we contributed specific advice or analysis for your work.
          </p>
        </div>
      </main>
    </Layout>
  )
}

export default Downloads
