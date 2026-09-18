# AWS and container image versions

The CI and ECR publishing workflows pull the amd64 OpenTOPAS base image from Docker Hub:
`opentopas/opentopas:v4.3.0-geant4-11.3.2-amd64`. Geant4 11.4.2 is currently incompatible with public TOPAS-nBio chemistry.
The Dockerfile receives this selection through the `OPENTOPAS_IMAGE` build argument.
The CI wrapper explicitly forwards it through the launcher's `-image=` option.

The default ECR TOPAS-nBio output image tag is `topas-v4.3.0-geant4-11.3.2-amd64`.
Here, `topas-v4.3.0` identifies the OpenTOPAS base version, not a TOPAS-nBio release number.
The checked-out TOPAS-nBio source is compiled into that image. For distinct nBio releases,
supply distinct output tags and update the Batch job definition to match.

- The AWS workflow pulls its OpenTOPAS base image from Docker Hub through `opentopas/opentopas`.
- AWS Batch template: `public.ecr.aws/q0u0d8d4/topas-nbio:topas-v4.3.0-geant4-11.3.2-amd64`
- Postprocessing: `public.ecr.aws/q0u0d8d4/topas-nbio:postprocessing`

The ECR publishing destination is configured by the GitHub secret `ECR_REPOSITORY_URI`.
It must match the repository in the Batch templates, or those templates must be
adjusted to the chosen destination. The secret should contain the repository URI without a tag.

Publish the OpenTOPAS base image first, then run `.github/workflows/build-push-ecr.yml`
to build and publish the combined TOPAS-nBio image to AWS ECR. No separate TOPAS-nBio
Docker Hub repository is required for CI or AWS Batch.
The ECR workflow remains amd64-only and also publishes a `latest` alias for existing users;
the simulation Batch template selects the explicit versioned tag.
Existing ECR images are retained, so publishing does not remove images used by older job definitions.

After publishing to ECR, register the updated `batch-job-definition.json` in AWS Batch.
Editing the local file alone does not change the registered AWS configuration.
`topas_submit.sh` submits to the named job definition; set `JOB_DEFINITION` to
`topas-nbio-job:<revision>` when a specific registered revision is required.
The postprocessing image and job definition do not depend on Geant4 and keep their existing tag.

## Optional Docker Hub publishing

`.github/workflows/build-push-dockerhub.yml` is a separate, manually triggered workflow.
It uses the same pinned OpenTOPAS base image with Geant4 11.3.2 and builds amd64 only.
Its intended destination is `opentopas/topas-nbio`, with the input tag defaulting to
`latest`; postprocessing uses `postprocessing`. This is a configured destination,
not a claim that the repository or images have already been published.

Docker Hub can create a missing repository on the first `docker push`, using the
namespace's default privacy settings. The configured `DOCKERHUB_USERNAME` and
`DOCKERHUB_TOKEN` must authorize publishing and, if needed, repository creation in
the `opentopas` organization. Organization restrictions can prevent automatic creation.
Alternatively, an organization administrator can create the repository beforehand
and grant the publishing credentials access. See [Docker Hub namespace settings](https://docs.docker.com/docker-hub/settings/).

This optional workflow publishes the combined TOPAS-nBio image to its own repository;
it does not replace the base image in `opentopas/opentopas` or change AWS Batch's ECR selection.
