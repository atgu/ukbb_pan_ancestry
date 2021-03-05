import { makeStyles, Switch, SwitchProps} from "@material-ui/core"
import React from "react"
import {  PopulationCode, populationColorMapping } from "./populations";

// Adapted from https://material-ui.com/components/switches/#customized-switches
// This lighten the switch's color by 50% when it's not selected:
const perPopulationStyles = {}
for (const [populationCode, color] of populationColorMapping.entries()) {
  perPopulationStyles[`${populationCode}SwitchBase`] = {
    '&$checked': {
      color: color,
    },
    '&$checked + $track': {
      backgroundColor: color,
    },
  }
}

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
