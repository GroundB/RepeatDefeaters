# Configs - Workflow specific configs

This folder contains configuration specific to
the operation of the workflow and it's tools.
Workflow analysis parameters are instead kept
in their own folders under the `analyses` folder.

The `compute_resources.config` file contains
the specification of cluster resources to request
from the cluster queue manager.

The `modules.config` file contains default
parameters for tools used in the workflow.
They can be overridden in a custom config to
provide alternative parameters to specific tools.
