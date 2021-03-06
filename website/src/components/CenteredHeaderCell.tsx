import { makeStyles } from "@material-ui/core"
import React from "react"
import { Tooltip } from "@material-ui/core";

const useStyles = makeStyles({
  root: {
    display: "flex",
    justifyContent: "center",
    alignItems: "center",
    height: "100%",
  }
} )
interface Props {
  column: {columnTitle: string}
}

export const CenteredHeaderCell = ({column: {columnTitle}}: Props): React.ReactNode => {
  const classes = useStyles()
  return (
    <div className={classes.root}><div>{columnTitle}</div></div>
  )
}

export const SaigeHeritabilityHeaderCell = ({column: {columnTitle}}: Props) => {
  const classes = useStyles()
  const tooltipText = "Warning: heritability computed by SAIGE is downwardly biased. See About > Technical Details > QC for details"
  return (
    <Tooltip title={tooltipText} placement="top">
      <div className={classes.root}><div>{columnTitle}</div></div>
    </Tooltip>
  )
}
