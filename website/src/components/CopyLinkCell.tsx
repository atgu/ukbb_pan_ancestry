import React from "react"
import { IconButton, } from '@material-ui/core';
import {  FileCopy as FileCopyIcon } from "@material-ui/icons";

interface Props {
  value: string
}

export const CopyLinkCell = ({value}: Props) => {
  const onClick = () => {
    // TODO: actually copy to clipboard here
  }
  return (
    <div style={{display: "flex", justifyContent: "center", alignItems: "center"}}>
      <IconButton color="primary" onClick={onClick}>
        <FileCopyIcon/>
      </IconButton>
    </div>
  )
}

export const copyLinkCellMaxWidth = 65
