/* a few ASCII definitions */

#define NULL_CHAR		  0
#define END_OF_STRING		  0
#define BELL			  7	/* ^G */
#define BACKSPACE		  8	/* ^H */
#define TAB			  9	/* ^I */
#define LINE_FEED		 10	/* ^J */
#define VERTICAL_TAB		 11	/* ^K */
#define FORM_FEED		 12	/* ^L */
#define CARRIAGE_RETURN		 13	/* ^M */
#define ESC			 27	/* ^[ */
#define SPACE			 32
#define END_OF_FILE		255

/* Declare some VT100 escape sequence strings; define with init_esc_seq(). */

char CLEAR_SCREEN[10], CLEAR_EOS[10], CLEAR_EOL[10], 
	CURSOR_SAVE[10], CURSOR_RESTORE[10], CURSOR_HOME[10], 
	CURSOR_UP[10], CURSOR_DOWN[10], CURSOR_RIGHT[10], CURSOR_LEFT[10],
	LINE_INSERT[10], 
	CURSOR_OFF[10], CURSOR_ON[10], 
	GRAPHICS_ON[10], GRAPHICS_OFF[10],
	TEXT_NORMAL[10], TEXT_BOLD[10], TEXT_UNDERLINE[10], TEXT_BLINK[10], 
	TEXT_REVERSE[10], TEXT_BOLD_OFF[10], TEXT_UNDERLINE_OFF[10], 
	TEXT_BLINK_OFF[10], TEXT_REVERSE_OFF[10],
	FONT_SELECT_DHT[10], FONT_SELECT_DHB[10], 
	FONT_SELECT_SW[10], FONT_SELECT_DW[10],
	SCROLL_RESTORE[10];

/*****************************************************************************
 **									    **
 **			     I N I T _ E S C _ S E Q			    **
 **									    **
 **   Initialize VT100 escape sequence strings for use in screen output.    **
 **									    **
 *****************************************************************************/
void init_esc_seq ()
{
  /* all of these are global character arrays declared in CHARDEFS.H */
  sprintf (CLEAR_SCREEN,    "%c[2J", ESC);  /* erase entire screen */
  sprintf (CLEAR_EOS,       "%c[J",  ESC);  /* erase to end of screen */
  sprintf (CLEAR_EOL,       "%c[K",  ESC);  /* erase to end of line */
  sprintf (CURSOR_SAVE,     "%c[s",  ESC);  /* save cursor position */
  sprintf (CURSOR_RESTORE,  "%c[u",  ESC);  /* restore saved cursor position */
  sprintf (CURSOR_HOME,     "%c[H",  ESC);  /* ESC [ row ; col H is general */
  sprintf (CURSOR_UP,       "%c[1A", ESC);  /* ESC [ n A moves n lines up */
  sprintf (CURSOR_DOWN,     "%c[1B", ESC);  /* ESC [ n B moves n lines down */
  sprintf (CURSOR_RIGHT,    "%c[1C", ESC);  /* ESC [ n C moves n lines right */
  sprintf (CURSOR_LEFT,     "%c[1D", ESC);  /* ESC [ n D moves n lines left */
  sprintf (LINE_INSERT,     "%c[1L", ESC);  /* push down lines below this one*/
  sprintf (CURSOR_OFF,      "%c[?25l", ESC);
  sprintf (CURSOR_ON,       "%c[?25h", ESC);
  sprintf (GRAPHICS_ON,     "%c(0",  ESC);  /* character graphics switch */
  sprintf (GRAPHICS_OFF,    "%c(B",  ESC);
  sprintf (TEXT_NORMAL,     "%c[0m", ESC);  /* resets all text modes */
  sprintf (TEXT_BOLD,       "%c[1m", ESC);
  sprintf (TEXT_UNDERLINE,  "%c[4m", ESC);
  sprintf (TEXT_BLINK,      "%c[5m", ESC);
  sprintf (TEXT_REVERSE,    "%c[7m", ESC);
  sprintf (TEXT_BOLD_OFF,      "%c[22m", ESC);  /* turn off specific modes */
  sprintf (TEXT_UNDERLINE_OFF, "%c[24m", ESC);
  sprintf (TEXT_BLINK_OFF,     "%c[25m", ESC);
  sprintf (TEXT_REVERSE_OFF,   "%c[27m", ESC);
  sprintf (FONT_SELECT_DHT, "%c#3",  ESC); /* double height+width top line */
  sprintf (FONT_SELECT_DHB, "%c#4",  ESC); /* double height+width bottom line*/
  sprintf (FONT_SELECT_SW,  "%c#5",  ESC); /* single height+width top line */
  sprintf (FONT_SELECT_DW,  "%c#6",  ESC); /* double width, single height ln */
  sprintf (SCROLL_RESTORE,  "%c[1;24r%c[24;1H", ESC, ESC); /* set normal scrl*/

		/* the general command to set a scrolling window and move the
		   cursor to the bottom line in that scrolling window is:
			ESC [ top ; bot r ESC [ bot ; 1 H */

  /*
     some online references:
	(obtained by searching google.com on "vt100 escape sequences")

	http://www.csd.uch.gr/~athaks/manuals/vt100esc.html
	http://www.termsys.demon.co.uk/vtansi.htm
	http://www.coe.uncc.edu/~danderse/www/vt100-ctrlcode.html
	http://vt100.net/terminals_faq
  */
}
