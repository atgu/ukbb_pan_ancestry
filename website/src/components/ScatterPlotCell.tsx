import React, { useRef, useState } from "react"
import {commonPopulations, PerPopulationMetrics, populationColorMapping} from "./populations";
import {  scaleLinear } from "d3-scale";
import { makeStyles, Tooltip } from "@material-ui/core";
import { ColumnInstance } from "react-table";
import { Datum } from "./types";
import {  height, width } from "./BarChartCell";

type XAxisConfig = {
  isShown: false
} | {
  isShown: true,
  yIntercept: number
}
interface Props {
  value: PerPopulationMetrics
  column: ColumnInstance<Datum> & {
    minValue: number,
    maxValue: number,
    xAxisConfig: XAxisConfig
    labelFormatter: (value: number) => string
  }
}

const displayedPopulations = commonPopulations

const margins = {
  left: 10,
  right: 10,
  top: 13,
  bottom: 10,
}
const dotDiameter = 8
const horizontalDistanceBetweenDots = (width - margins.left - margins.right) / commonPopulations.length
const xScale = scaleLinear<number, number>().domain([
  0, commonPopulations.length - 1
]).range([
  margins.left + horizontalDistanceBetweenDots / 2, width - margins.right - horizontalDistanceBetweenDots / 2
])

const useStyles = makeStyles({
  dotContainer: {
    position: "absolute",
    top: "0",
  },
  barLabel: {
    position: "absolute",
    fontSize: "0.7rem",
    transform: "translate(-50%, -50%)",
  },
})

export const ScatterPlotCell = ({value, column}: Props) => {
  const classes = useStyles()
  const {minValue, maxValue, xAxisConfig, labelFormatter} = column
  const yScale = scaleLinear<number, number>().
    domain([minValue, maxValue]).
    range([height - margins.bottom, margins.top]).
    clamp(true)
  const bars = displayedPopulations.map((population, index) => {
    if (population === undefined) {
      return null
    } else {
      const numericValue = value.get(population)!
      if (numericValue === undefined) {
        return null
      } else {
        const color = populationColorMapping.get(population)
        const topYCoord = yScale(numericValue)
        const label = (
          <div
            className={classes.barLabel}
            style={{ top: `${topYCoord - 9}px`, }}
          >{labelFormatter(numericValue)}</div>
        )
        return (
          <div className={classes.dotContainer} key={population} style={{
            left: `${xScale(index)}px`,
            height: `${height}px`,
          }}>
            <div key={population}
              style={{
                position: "absolute",
                width: `${dotDiameter}px`,
                height: `${dotDiameter}px`,
                left: "0",
                top: `${topYCoord}px`,
                transform: `translate(-50%, -50%)`,
                borderRadius: "50%",
                backgroundColor: color,
              }}
            >
            </div>
            {label}
          </div>
        )
      }
    }
  })
  const xAxisElem = (xAxisConfig.isShown === true) ? (
    <div style={{
      position: "absolute",
      left: "0px",
      right: "0px",
      height: "1px",
      top: `${yScale(xAxisConfig.yIntercept)}px`,
      backgroundColor: "black",
    }} />
  ) : null

  return (
    <div style={{width: `${width}px`, height: `${height}px`, position: "relative"}}>
      {xAxisElem}
      {bars}
    </div>
  )
}
