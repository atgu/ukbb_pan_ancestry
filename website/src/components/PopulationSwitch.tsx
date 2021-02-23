import { makeStyles, Switch, SwitchProps} from "@material-ui/core"
import React from "react"
import {  lighten } from "polished";
import {  PopulationCode, populationColorMapping } from "./populations";


// Adapted from https://material-ui.com/components/switches/#customized-switches
// This lighten the switch's color by 50% when it's not selected:
const perPopulationStylePairs = [...populationColorMapping.entries()].map(([populationCode, color]) => (
  [`${populationCode}SwitchBase`, {
    color: lighten(0.2, color),
    '&$checked': {
      color: color,
    },
    '&$checked + $track': {
      backgroundColor: color,
    },
  }]
))
const perPopulationStyles = Object.fromEntries(perPopulationStylePairs)

const useStyles = makeStyles({
  ...perPopulationStyles,
  checked: {},
  track: {}
})

interface Props {
  population: PopulationCode
}
export const PopulationSwitch = (props: Props & SwitchProps) => {
  const classes = useStyles()
  const {population, ...rest} = props;
  return (
    <Switch
      classes={{
        switchBase: classes[`${population}SwitchBase`],
        checked: classes[`checked`],
        track: classes[`track`],
      }}
      {...rest}
    />
  )
}
