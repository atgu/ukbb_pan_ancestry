import React from 'react'
import Layout from '@theme/Layout'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import {  PhenotypesPageContent } from "../components/PhenotypesPageContent";


const Phenotypes = () => {
  const context = useDocusaurusContext()
  const { siteConfig = {} } = context
  return (
    <Layout title={`${siteConfig.title}`} description="Phenotypes">
      <PhenotypesPageContent/>
    </Layout>
  )
}

export default Phenotypes
