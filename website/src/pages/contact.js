import React from 'react'
import classnames from 'classnames'
import Layout from '@theme/Layout'
import Link from '@docusaurus/Link'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import useBaseUrl from '@docusaurus/useBaseUrl'
import styles from './styles.module.css'

function Contact() {
  const context = useDocusaurusContext()
  const { siteConfig = {} } = context
  return (
    <Layout title={`${siteConfig.title}`} description="Contact">
      <header className={classnames('hero', styles.heroBanner)}>
        <div className="container">
          <h1>Contact</h1>
        </div>
      </header>
      <main>
        <div className="container">
          <p>
            For questions, data problems, or feature suggestions, <a href="mailto:ukb.diverse.gwas@gmail.com">email us</a>.
          </p>
          <p>
            For details about the code to run this analysis, see our <a href="https://github.com/atgu/ukbb_pan_ancestry">Github</a>.
          </p>
        </div>
      </main>
    </Layout>
  )
}

export default Contact
