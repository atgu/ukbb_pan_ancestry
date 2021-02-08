import { makeStyles, Tooltip } from "@material-ui/core"
import React from "react"

const useStyles = makeStyles(() => ({
  cell: {
    whiteSpace: "nowrap",
    overflow: "hidden",
    textOverflow: "ellipsis",
  }
}))

interface Props {
  value: string
}

export const TruncatedTextCell = ({value}: Props) => {
  const classes = useStyles()
  return (
    <Tooltip title={value} placement="top-start">
      <div className={classes.cell}>
        {value}
      </div>
    </Tooltip>
  )
}

export const downloadLinkCellMaxWidth = 50
