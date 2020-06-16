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
      <header className={classnames('hero', styles.heroBanner)}>
        <div className="container">
          <h1>Downloads</h1>
        </div>
      </header>
      <main>
        <div className="container">
          <p class="">
            The dataset is available for download in two main formats:
          </p>
          <ul>
            <li>
              <b>Per-phenotype flat files</b>: for most analyses of one or a few phenotypes, we suggest using the
              per-phenotype flat files, available freely on Dropbox. More information on the file formats is
              available in the <Link to={useBaseUrl('docs/per-phenotype-files')}>Technical Details</Link>.
              The phenotype manifest (browse on <a href="https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=30994804" target="_blank">Google Sheets</a> or
              download on <a href="https://www.dropbox.com/s/18p4lj3finj11oh/phenotype_manifest.tsv.bgz?dl=0">Dropbox</a>) contains the location and
              detailed information of all phenotypes for which GWAS were run.
            </li>
            <li>
              <b>Hail format</b>: For large-scale analyses of many phenotypes, we provide the full dataset
              in <Link to={useBaseUrl('docs/hail-format')}>Hail MatrixTable format</Link> on Google Cloud.
            </li>
          </ul>
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
