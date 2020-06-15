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
          <p class="hero__description text--center">
            Email us: <a href="mailto:ukb.diverse.gwas@gmail.com">ukb.diverse.gwas@gmail.com</a>
          </p>
        </div>
      </main>
    </Layout>
  )
}

export default Contact
