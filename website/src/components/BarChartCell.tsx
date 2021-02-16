import React, { useRef, useState } from "react"
import {commonPopulations, PerPopulationMetrics, populationColorMapping} from "./populations";
import {  scaleLinear } from "d3-scale";
import { makeStyles, Tooltip } from "@material-ui/core";
import { ColumnInstance } from "react-table";
import { Datum } from "./types";


export const height = 60
export const width = 200;

interface Props {
  value: PerPopulationMetrics
  column: ColumnInstance<Datum> & {minValue: number, maxValue: number}
}

const displayedPopulations = commonPopulations

const margins = {
  left: 10,
  right: 10,
  top: 5,
  bottom: 5,
}
const barWidth = 20
const horizontalDistanceBetweenBars = (width - margins.left - margins.right  - commonPopulations.length * barWidth) / commonPopulations.length
const xScale = scaleLinear<number, number>().domain([
  0, commonPopulations.length - 1
]).range([
  margins.left + horizontalDistanceBetweenBars / 2, width - margins.right - horizontalDistanceBetweenBars / 2
])

export const BarChartCell = ({value, column}: Props) => {
  const {minValue, maxValue} = column
  const yScale = scaleLinear<number, number>().domain([minValue, maxValue]).range([height - margins.bottom, margins.top])
  const bars = displayedPopulations.map((population, index) => {
    if (population === undefined) {
      return null
    } else {
      const numericValue = value.get(population)
      const color = populationColorMapping.get(population)
      let barHeight: number, y: number
      if (numericValue === 0) {
        const epsilonHeight = 1
        barHeight = epsilonHeight
        y = height - margins.bottom - epsilonHeight
      } else {
        const topYCoord = yScale(numericValue)
        barHeight = yScale(minValue) - topYCoord
        y = topYCoord
      }

      return (
        <Tooltip title={`${population}: ${numericValue}`} placement="top">
          <div key={population}
            style={{
              position: "absolute",
              width: `${barWidth}px`,
              height: `${barHeight}px`,
              bottom: `${margins.bottom}px`,
              left: `${xScale(index) - barWidth / 2}px`,
              backgroundColor: color,
            }}
          />
        </Tooltip>
      )
    }
  })

  return (
    <div style={{width: `${width}px`, height: `${height}px`, position: "relative"}}>
      {bars}
    </div>
  )
}
