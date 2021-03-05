import React, { useEffect, useRef } from "react"
import Mousetrap from "mousetrap"
import { BaseTextFilter } from "./BaseTextFilter"


export const GlobalFilter = ({preGlobalFilteredRows, globalFilter, setGlobalFillter}) => {
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
  const count = preGlobalFilteredRows.length
  return (
    <BaseTextFilter
      ref={inputRef}
      filterValue={globalFilter}
      setFilterValue={setGlobalFillter}
      label={`Search all fields in ${count} records`}
    />
  )
}
