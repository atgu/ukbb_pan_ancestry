import { makeStyles } from "@material-ui/core";
import { FiberManualRecord } from "@material-ui/icons";
import React from "react"
import {  populationColorMapping, commonPopulations as expectedPopulations } from "./populations";

const useStyles = makeStyles({
  grid: {
    display: "grid",
    gridTemplateColumns: "repeat(3, 1fr)",
    columnGap: "5px",
    gridTemplateRows: "repeat(2, 15px)",
    justifyContent: "center",
    alignItems: "center",
  },
  population: {
    display: "flex",
    alignItems: "center",
  },
  icon: {
    fontSize: "1rem",
  },
  populationName: {
    fontSize: "0.75rem",
  }
})
export const PopuplationHeader = ({column}) => {
  const classes = useStyles()
  const populationElems = expectedPopulations.map(pop => {
    return (
      <div key={pop} className={classes.population}>
        <FiberManualRecord  style={{color: populationColorMapping.get(pop)}} className={classes.icon}/>
        <div className={classes.populationName}>{pop}</div>
      </div>
    )
  })
  return (
    <div {...column.getHeaderProps()}>
      Populations
      <div className={classes.grid}>
        {populationElems}
      </div>
    </div>
  )
}
