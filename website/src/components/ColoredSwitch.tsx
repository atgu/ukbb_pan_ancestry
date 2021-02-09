import { ButtonProps, makeStyles, Switch, SwitchProps, withStyles } from "@material-ui/core"
import React from "react"

export const NewColoredSwitch = withStyles({
  switchBase: (props: Props) => ({
    color: props.color,
    '&$checked': {
      color: props.color,
    },
    '&$checked + $track': {
      backgroundColor: props.color,
    },
  }),
  checked: {},
  track: {},
})(Switch);

interface Props {
  color: string
}
