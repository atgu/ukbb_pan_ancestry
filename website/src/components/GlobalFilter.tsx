import React, { useEffect, useRef, useState } from "react"
import { makeStyles, TextField } from "@material-ui/core"
import {  useAsyncDebounce } from "react-table";
import Mousetrap from "mousetrap"

const useStyles = makeStyles(() => ({
  labelUnfocused: {
    whiteSpace: "nowrap",
    overflow: "hidden",
    textOverflow: "ellipsis",
    maxWidth: "100%"
  },
  labelFocused: {
    whiteSpace: "unset",
    overflow: "unset",
    textOverflow: "unset",
    maxWidth: "unset",
  }
}))

export const GlobalFilter = ({preGlobalFilteredRows, globalFilter, setGlobalFillter}) => {
  const classes = useStyles()
  const count = preGlobalFilteredRows.length
  const [value, setValue] = useState(globalFilter)
  const onChange = useAsyncDebounce(val => {
    setGlobalFillter(value || "")
  }, 250)
  const inputRef = useRef<HTMLInputElement | null>(null)
  const keyboardHandlerRef = useRef((e: KeyboardEvent) => {
    const {current: inputElem} = inputRef;
    if (inputElem !== null) {
      // Need to avoid triggering the native browser search box:
      e.stopPropagation()
      e.preventDefault()
      inputElem.focus()
    }
  })
  useEffect(() => {
    // "mod+f" means "cmd+f" on Mac and "ctrl+f" on Windows:
    const shortcuts = ["mod+f"]

    Mousetrap.bind(shortcuts, keyboardHandlerRef.current)
    return () => {
      Mousetrap.unbind(shortcuts)
    }
  }, [])
  return (
    <TextField
      type="search"
      fullWidth={true}
      value={value}
      onChange={e => {
        setValue(e.target.value)
        onChange(e.target.value)
      }}
      inputRef={inputRef}
      InputLabelProps={{
        classes: {
          root: classes.labelUnfocused,
          focused: classes.labelFocused,
        }
      }}
      label={ `Search all fields in ${count} records` }
    />

  )
}
