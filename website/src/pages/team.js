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
          <p>
            Analysis team:
            <ul>
              <li>Konrad Karczewski</li>
              <li>Elizabeth Atkinson</li>
              <li>Nikolas Baya</li>
              <li>Masahiro Kanai</li>
              <li>Alicia Martin</li>
              <li>Patrick Turley</li>
              <li>Gopal Sarma</li>
              <li>Raymond Walters</li>
              <li>Duncan Palmer</li>
              <li>Nathan Cheng</li>
              <li>Wei Zhou</li>
              <li>Rahul Gupta</li>
            </ul>
            Ethics team
            <ul>
              <li>Patrick Turley</li>
              <li>Elizabeth Atkinson</li>
              <li>Shawneequa Callier</li>
              <li>Alicia Martin</li>
            </ul>
            Website team:
            <ul>
              <li>Matthew Solomonson</li>
              <li>Konrad Karczewski</li>
            </ul>
            Data coordinators
            <ul>
              <li>Sam Bryant</li>
              <li>Claire Churchhouse</li>
              <li>Caroline Cusick</li>
            </ul>
            Hail
            <ul>
              <li>Jacqueline Goldstein</li>
              <li>Daniel King</li>
              <li>Cotton Seed</li>
            </ul>
            PIs
            <ul>
              <li>Alicia Martin</li>
              <li>Hilary Finucane</li>
              <li>Benjamin Neale</li>
              <li>Mark Daly</li>
            </ul>
          </p>
        </div>
      </main>
    </Layout>
  )
}

export default Contributors
