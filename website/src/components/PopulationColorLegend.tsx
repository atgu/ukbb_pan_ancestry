
import { makeStyles } from "@material-ui/core";
import { FiberManualRecord } from "@material-ui/icons";
import React from "react"
import {  populationColorMapping, commonPopulations as expectedPopulations } from "./populations";

const useStyles = makeStyles({
  grid: {
    display: "grid",
    gridTemplateColumns: "repeat(3, 1fr)",
    gridTemplateRows: "repeat(2, 1fr)",
    justifyContent: "center",
    alignItems: "center",
    alignContent: "center",
    flexGrow: 1,
  },
  population: {
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
  },
})
interface Props {
  labelFontSize: string | undefined
  iconFontSize: string | undefined
  gridJustifyItems: string | undefined
}
export const PopulationColorLegend = ({labelFontSize, iconFontSize, gridJustifyItems}: Props) => {
  const classes = useStyles()
  const populationElems = expectedPopulations.map(pop => {
    return (
      <div key={pop} className={classes.population}>
        <FiberManualRecord  style={{
            color: populationColorMapping.get(pop),
            fontSize: iconFontSize,
        }}/>
        <div style={{ fontSize: labelFontSize}}>{pop}</div>
      </div>
    )
  })
  return (
    <div className={classes.grid} style={{justifyItems: gridJustifyItems}}>
      {populationElems}
    </div>
  )

}
