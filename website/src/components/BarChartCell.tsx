import React  from "react"
import {commonPopulations, PerPopulationMetrics, populationColorMapping} from "./populations";
import {  scaleLinear } from "d3-scale";
import { makeStyles} from "@material-ui/core";
import { ColumnInstance } from "react-table";
import { Datum } from "./types";

export const height = 80
export const width = 200;



interface Props {
  value: PerPopulationMetrics
  column: ColumnInstance<Datum> & {
    labelFormatter: (value: number) => string
    // Values greater than this threshold are shown with a broken bar:
    upperThreshold: number
  }
}

const displayedPopulations = commonPopulations

const margins = {
  left: 10,
  right: 10,
  top: 13,
  bottom: 5,
}
// Bars representing values above the `upperThreshold` will extend this much above the
const thresholdExcessVerticalSpace = 10;

const barWidth = 20
const horizontalDistanceBetweenBars = (width - margins.left - margins.right  - commonPopulations.length * barWidth) / commonPopulations.length
const xScale = scaleLinear<number, number>().domain([
  0, commonPopulations.length - 1
]).range([
  margins.left + horizontalDistanceBetweenBars / 2, width - margins.right - horizontalDistanceBetweenBars / 2
])

const useStyles = makeStyles({
  root: {
    position: "relative",
    width: `${width}px`,
    height: `${height}px`
  },
  brokenBarStripe: {
    position: "absolute",
    top: "15px",
    height: "3px",
    left: "0",
    right: "0",
    backgroundColor: "white",
    transform: "skewY(-20deg)",
  },
  barLabel: {
    position: "absolute",
    fontSize: "0.7rem",
    transform: "translate(-50%, -50%)",
  },
  barContainer: {
    position: "absolute",
    top: "0",
  },
  bar: {
    position: "absolute",
    left: 0,
    transform: "translateX(-50%)",
    width: `${barWidth}px`,
    bottom: `${margins.bottom}px`,
  }
})

export const BarChartCell = ({value, column}: Props): React.ReactNode => {
  const classes = useStyles()
  const {upperThreshold, labelFormatter} = column
  const minValue = 0
  const yScale = scaleLinear<number, number>().
    domain([minValue, upperThreshold]).
    range([height - margins.bottom, margins.top + thresholdExcessVerticalSpace])
  const bars = displayedPopulations.map((population, index) => {
    if (population === undefined) {
      return null
    } else {
      const numericValue = value.get(population)
      if (numericValue === undefined) {
        return null
      } else {
        const color = populationColorMapping.get(population)
        let barHeight: number, labelTop: number, diagonalStripe: React.ReactNode
        if (numericValue === 0) {
          const epsilonHeight = 1
          barHeight = epsilonHeight
          diagonalStripe = null
          labelTop = height - margins.bottom - barHeight
        } else if (numericValue > upperThreshold) {
          barHeight = yScale(minValue) - margins.top
          labelTop = margins.top
          diagonalStripe = (
            <div className={classes.brokenBarStripe} />
          )
        } else {
          const topYCoord = yScale(numericValue)
          barHeight = yScale(minValue) - topYCoord
          labelTop = topYCoord
          diagonalStripe = null
        }
        const labelText =  (numericValue === 0) ? 0 : labelFormatter(numericValue)
        const label = (
          <div
            className={classes.barLabel}
            style={{ top: `${labelTop - 6}px`, }}
          >{labelText}</div>
        )

        return (
          <div className={classes.barContainer} key={population} style={{
            left: `${xScale(index)}px`,
            height: `${height}px`,
          }}>
            <div
              className={classes.bar}
              style={{
                height: `${barHeight}px`,
                backgroundColor: color,
              }}
            >
              {diagonalStripe}
            </div>
            {label}
          </div>
        )
      }
    }
  })

  return (
    <div className={classes.root}>
      {bars}
    </div>
  )
}
