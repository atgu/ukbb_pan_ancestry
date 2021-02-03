import React from "react"
import {  populationColorMapping, commonPopulations as expectedPopulations } from "./populations";
import styles from "../pages/phenotypes.module.css"
import { FiberManualRecord } from "@material-ui/icons";
import { Tooltip } from "@material-ui/core"

export const PopulationCell = ({value}) => {
  const actualPopulations = value ? new Set(value.split(",")) : new Set()
  const populationElems = expectedPopulations.map(pop => {
    if (actualPopulations.has(pop) ) {
      return (
        <Tooltip title={pop} placement="top">
          <FiberManualRecord style={{color: populationColorMapping.get(pop)}}/>
        </Tooltip>
      )
    } else {
      return <div/>
    }
  })

  return (
    <div className={styles.populationCell}>
      {populationElems}
    </div>
  )
}
