import React from "react"
import {  populationColorMapping, commonPopulations as expectedPopulations } from "./populations";
import { FiberManualRecord } from "@material-ui/icons";
import { makeStyles, Tooltip } from "@material-ui/core"

const useStyles = makeStyles({
  root: {
    display: "grid",
    gridTemplateColumns: "repeat(3, 25px)",
    gridTemplateRows: "repeat(2, 25px)",
    justifyContent: "center",
    alignItems: "center",
  }
})

export const PopulationCell = ({value}) => {
  const classes = useStyles()
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
    <div className={classes.root}>
      {populationElems}
    </div>
  )
}
