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
      }),
    [isDarkTheme],
  );
  return (
    <ThemeProvider theme={materialUiTheme}>
      <PhenotypesPageContent/>
    </ThemeProvider>
  )
}
