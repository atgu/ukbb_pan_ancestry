import React from "react"
import { createMuiTheme, ThemeProvider } from "@material-ui/core"
import useThemeContext from '@theme/hooks/useThemeContext';
import { PhenotypesPageContent } from "./PhenotypesPageContent";


export const PhenotypesPageContentWrapper = () => {
  const {isDarkTheme} = useThemeContext();

  const materialUiTheme = React.useMemo(
    () =>
      createMuiTheme({
        palette: {
          type: isDarkTheme ? 'dark' : 'light',
        },
        overrides: {
          MuiTableCell: {
            root: {
              padding: "0",
            }
          }
        }
      }),
    [isDarkTheme],
  );
  return (
    <ThemeProvider theme={materialUiTheme}>
      <PhenotypesPageContent/>
    </ThemeProvider>
  )
}
