import React, { useRef, useState } from "react"
import {commonPopulations, PerPopulationMetrics, populationColorMapping} from "./populations";
import {  scaleLinear } from "d3-scale";
import { makeStyles, Tooltip } from "@material-ui/core";


export const height = 60
export const width = 200;

interface Props {
  value: PerPopulationMetrics
}

const minValue = 0
const maxValue = 1
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
const yScale = scaleLinear<number, number>().domain([minValue, maxValue]).range([height - margins.bottom, margins.top])

export const SaigeHeritabilityCell = ({value}: Props) => {
  const bars = displayedPopulations.map((population, index) => {
    const numericValue = value.get(population)
    if (population === undefined) {
      return null
    } else {
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
  const yCoordOfXAxis = yScale(minValue)
  return (
    <div style={{width: `${width}px`, height: `${height}px`, position: "relative"}}>
      {bars}
    </div>
  )
}
