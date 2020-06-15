import React from 'react'
import classnames from 'classnames'
import Layout from '@theme/Layout'
import Link from '@docusaurus/Link'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import useBaseUrl from '@docusaurus/useBaseUrl'
import styles from './styles.module.css'

function Contributors() {
  const context = useDocusaurusContext()
  const { siteConfig = {} } = context
  return (
    <Layout title={`${siteConfig.title}`} description="Contributors">
      <header className={classnames('hero', styles.heroBanner)}>
        <div className="container">
          <h1>Contributors</h1>
        </div>
      </header>
      <main>
        <div className="container">
          <p class="hero__description text--center">
            Konrad J. Karczewski, Elizabeth G. Atkinson, Masahiro Kanai, Nikolas Baya, Patrick Turley,
            Shawneequa Callier, Gopal Sarma, Raymond K. Walters, Duncan S. Palmer, Matthew Solomonson, Nathan
            Cheng, Rahul Gupta, Sam Bryant, Claire Churchhouse, Caroline Cusick, Jacqueline I. Goldstein,
            Daniel King, Wei Zhou, Cotton Seed, Mark J. Daly, Benjamin M. Neale, Hilary Finucane, Alicia R.
            Martin.
          </p>
        </div>
      </main>
    </Layout>
  )
}

export default Contributors
