import { makeStyles, Tooltip } from "@material-ui/core"
import React from "react"

const useStyles = makeStyles(() => ({
  cell: {
    // https://dropshado.ws/post/1015351370/webkit-line-clamp
    overflow: "hidden",
    textOverflow: "ellipsis",
    display: "-webkit-box",
    "-webkit-line-clamp": 2,
    "-webkit-box-orient": "vertical",
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
