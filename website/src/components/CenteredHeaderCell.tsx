import { makeStyles } from "@material-ui/core"
import React from "react"

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

export const CenteredHeaderCell = ({column: {columnTitle}}: Props) => {
  const classes = useStyles()
  return (
    <div className={classes.root}><div>{columnTitle}</div></div>
  )
}
