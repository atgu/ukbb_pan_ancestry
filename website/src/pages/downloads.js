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
          <p class="hero__description text--center">
            The downloads go here.
          </p>
        </div>
      </main>
    </Layout>
  )
}

export default Downloads
