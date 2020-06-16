import React from 'react'
import classnames from 'classnames'
import Layout from '@theme/Layout'
import Link from '@docusaurus/Link'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import useBaseUrl from '@docusaurus/useBaseUrl'
import styles from './styles.module.css'

function Results() {
  const context = useDocusaurusContext()
  const { siteConfig = {} } = context
  return (
    <Layout title={`${siteConfig.title}`} description="Results">
      <header>
        <div className="container">
          <h1 className="page-title">Results</h1>
        </div>
      </header>
      <main>
        <div className="container">
          <p className="hero__description text--center">
            Coming soon.
          </p>
        </div>
      </main>
    </Layout>
  )
}

export default Results
