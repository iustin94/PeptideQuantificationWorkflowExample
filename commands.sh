
#Script to run x!tandem search engine on spectrum files and generate peptishaker project and report

set -e #Exit if any command fails"

SEARCHGUI_BIN=$1
PEPTIDESHAKER_BIN=$2
SPECTRUM_FILES=$3
PARAMETERS_FILE=$4

BASE_DIR=$(dirname "$0")

#Make folders in current execution folder
mkdir -p "SearchGUI"
mkdir -p "SearchGUI/Results"
mkdir -p "SearchGUI/Logs"

mkdir -p "PeptideShaker"
mkdir -p "PeptideShaker/Results"
mkdir -p "PeptideShaker/Logs"

#File paths parameters
SRC_GUI_RESULTS="$BASE_DIR/SearchGUI/Results"
SRC_GUI_LOGS="$BASE_DIR/SearchGUI/Logs"
PEP_SHAKE_RESULTS="$BASE_DIR/PeptideShaker/Results"
PEP_SHAKE_LOGS="$BASE_DIR/PeptideShaker/Logs"

echo "Running SearchGui ... "

java -cp "$SEARCHGUI_BIN/SearchGui-3.2.20.jar" "$SEARCHGUI_BIN/eu.isas.searchgui.cmd.SearchCLI" \
-spectrum_files $SPECTRUM_FILES \
-output_folder "$SRC_GUI_RESULTS" \
-id_params $PARAMETERS_FILE \
-xtandem "1" \
-decoy \
-log "$SRC_GUI_LOGS"

echo "Running PeptideShaker"

java -cp "$PEPTIDESHAKER_BIN/PeptideShaker-1.16.15.jar" "$PEPTIDESHAKER_BIN/eu.isas.peptideshaker.cmd.PeptideShakerCLI" \
-experiment "Malarya" \
-sample "Natalia" \
-replicate "1" \
-identification_files "$SRC_GUI_RESULTS/searchgui_out.zip"\
-out "$PEP_SHAKE_RESULTS/output.cpsx"\
-log "$PEP_SHAKE_LOGS"

echo "Generating PeptideShaker report" 

java -cp "$PEPTIDESHAKER_BIN/PeptideShaker-1.16.15.jar" "$PEPTIDESHAKER_BIN/eu.isas.peptideshaker.cmd.ReportCLI" \
-in "$PEP_SHAKE_RESULTS/output.cpsx"\
-out_reports "$PEP_SHAKE_RESULTS"\
-reports "3"\
-documentation "3"

