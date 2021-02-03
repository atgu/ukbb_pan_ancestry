import { Tooltip } from "@material-ui/core"
import React from "react"
import styles from "../pages/phenotypes.module.css"

interface Props {
  value: string
}

export const TruncatedTextCell = ({value}: Props) => {
  return (
    <Tooltip title={value} placement="top-start">
      <div className={styles.longTextCell}>
        {value}
      </div>
    </Tooltip>
  )
}

export const downloadLinkCellMaxWidth = 50
