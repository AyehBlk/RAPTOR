"""
Cloud Integration

Deploy and run RAPTOR pipelines on cloud platforms (AWS, Google Cloud, Azure).

Author: Ayeh Bolouki
Email: ayeh.bolouki@unamur.be
"""

import os
import json
from pathlib import Path
from typing import Dict, List, Optional
import logging
import subprocess

logger = logging.getLogger(__name__)


class CloudRunner:
    """
    Run RAPTOR pipelines on cloud platforms.
    
    Provides unified interface for deploying RNA-seq analysis pipelines
    on AWS, Google Cloud Platform, and Azure.
    
    Parameters
    ----------
    platform : str
        Cloud platform: 'aws', 'gcp', or 'azure'
    region : str, optional
        Cloud region (e.g., 'us-east-1' for AWS)
    credentials_file : str, optional
        Path to credentials file
    
    Examples
    --------
    >>> runner = CloudRunner(platform='aws', region='us-east-1')
    >>> runner.deploy_pipeline(
    ...     pipeline_id=3,
    ...     data_location='s3://my-bucket/data/',
    ...     output_location='s3://my-bucket/results/'
    ... )
    """
    
    def __init__(
        self,
        platform: str,
        region: Optional[str] = None,
        credentials_file: Optional[str] = None
    ):
        """Initialize cloud runner."""
        self.platform = platform.lower()
        self.region = region
        self.credentials_file = credentials_file
        
        # Validate platform
        if self.platform not in ['aws', 'gcp', 'azure']:
            raise ValueError(f"Unsupported platform: {platform}. Use 'aws', 'gcp', or 'azure'")
        
        # Set default regions
        if not self.region:
            defaults = {
                'aws': 'us-east-1',
                'gcp': 'us-central1',
                'azure': 'eastus'
            }
            self.region = defaults[self.platform]
        
        logger.info(f"CloudRunner initialized: {self.platform} ({self.region})")
        
        # Check credentials
        self._check_credentials()
    
    def _check_credentials(self):
        """Check if cloud credentials are configured."""
        if self.platform == 'aws':
            self._check_aws_credentials()
        elif self.platform == 'gcp':
            self._check_gcp_credentials()
        elif self.platform == 'azure':
            self._check_azure_credentials()
    
    def _check_aws_credentials(self):
        """Check AWS credentials."""
        # Check for AWS CLI
        try:
            result = subprocess.run(
                ['aws', '--version'],
                capture_output=True,
                text=True
            )
            logger.info(f"AWS CLI found: {result.stdout.strip()}")
        except FileNotFoundError:
            logger.warning("AWS CLI not found. Install from: https://aws.amazon.com/cli/")
            logger.warning("Pipeline deployment will not work without AWS CLI")
        
        # Check credentials
        if self.credentials_file and Path(self.credentials_file).exists():
            os.environ['AWS_SHARED_CREDENTIALS_FILE'] = self.credentials_file
            logger.info(f"Using credentials from: {self.credentials_file}")
        elif Path.home().joinpath('.aws', 'credentials').exists():
            logger.info("Using default AWS credentials")
        else:
            logger.warning("No AWS credentials found. Configure with: aws configure")
    
    def _check_gcp_credentials(self):
        """Check Google Cloud credentials."""
        # Check for gcloud CLI
        try:
            result = subprocess.run(
                ['gcloud', '--version'],
                capture_output=True,
                text=True
            )
            logger.info(f"gcloud CLI found: {result.stdout.split()[0]}")
        except FileNotFoundError:
            logger.warning("gcloud CLI not found. Install from: https://cloud.google.com/sdk")
        
        # Check credentials
        if self.credentials_file and Path(self.credentials_file).exists():
            os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = self.credentials_file
            logger.info(f"Using credentials from: {self.credentials_file}")
        elif os.getenv('GOOGLE_APPLICATION_CREDENTIALS'):
            logger.info("Using GCP credentials from environment")
        else:
            logger.warning("No GCP credentials found. Set GOOGLE_APPLICATION_CREDENTIALS")
    
    def _check_azure_credentials(self):
        """Check Azure credentials."""
        # Check for Azure CLI
        try:
            result = subprocess.run(
                ['az', '--version'],
                capture_output=True,
                text=True
            )
            logger.info("Azure CLI found")
        except FileNotFoundError:
            logger.warning("Azure CLI not found. Install from: https://docs.microsoft.com/cli/azure/install-azure-cli")
        
        # Azure credentials typically managed through az login
        logger.info("Azure credentials managed through: az login")
    
    def deploy_pipeline(
        self,
        pipeline_id: int,
        data_location: str,
        output_location: str,
        instance_type: Optional[str] = None,
        docker_image: Optional[str] = None,
        **kwargs
    ) -> Dict:
        """
        Deploy RAPTOR pipeline on cloud platform.
        
        Parameters
        ----------
        pipeline_id : int
            Pipeline ID (1-8)
        data_location : str
            Cloud storage location for input data
            (e.g., 's3://bucket/data/' for AWS)
        output_location : str
            Cloud storage location for outputs
        instance_type : str, optional
            Cloud instance type (e.g., 'c5.2xlarge' for AWS)
        docker_image : str, optional
            Docker image for RAPTOR (default: official image)
        **kwargs
            Additional platform-specific parameters
        
        Returns
        -------
        dict
            Deployment information including job ID and status
        """
        logger.info(f"Deploying Pipeline {pipeline_id} on {self.platform}")
        logger.info(f"Input: {data_location}")
        logger.info(f"Output: {output_location}")
        
        # Set default docker image
        if not docker_image:
            docker_image = "ayehblk/raptor:latest"
        
        # Platform-specific deployment
        if self.platform == 'aws':
            result = self._deploy_aws(
                pipeline_id, data_location, output_location,
                instance_type, docker_image, **kwargs
            )
        elif self.platform == 'gcp':
            result = self._deploy_gcp(
                pipeline_id, data_location, output_location,
                instance_type, docker_image, **kwargs
            )
        elif self.platform == 'azure':
            result = self._deploy_azure(
                pipeline_id, data_location, output_location,
                instance_type, docker_image, **kwargs
            )
        
        logger.info(f"Deployment successful: {result['job_id']}")
        
        return result
    
    def _deploy_aws(
        self,
        pipeline_id: int,
        data_location: str,
        output_location: str,
        instance_type: Optional[str],
        docker_image: str,
        **kwargs
    ) -> Dict:
        """Deploy on AWS using Batch or EC2."""
        logger.info("Deploying on AWS")
        
        # Default instance type
        if not instance_type:
            instance_type = "c5.2xlarge"  # 8 vCPUs, 16 GB RAM
        
        # Create job definition
        job_def = {
            'pipeline_id': pipeline_id,
            'platform': 'aws',
            'region': self.region,
            'instance_type': instance_type,
            'docker_image': docker_image,
            'data_location': data_location,
            'output_location': output_location
        }
        
        # For actual AWS Batch deployment (requires boto3)
        try:
            import boto3
            
            batch_client = boto3.client('batch', region_name=self.region)
            
            # Submit job to AWS Batch
            response = batch_client.submit_job(
                jobName=f'raptor-pipeline-{pipeline_id}',
                jobQueue=kwargs.get('job_queue', 'raptor-queue'),
                jobDefinition=kwargs.get('job_definition', 'raptor-job'),
                containerOverrides={
                    'command': [
                        'python', '/raptor/scripts/01_run_all_pipelines_python.py',
                        '--pipelines', str(pipeline_id),
                        '--data', data_location,
                        '--output', output_location
                    ],
                    'environment': [
                        {'name': 'RAPTOR_PIPELINE', 'value': str(pipeline_id)},
                        {'name': 'AWS_REGION', 'value': self.region}
                    ]
                }
            )
            
            job_id = response['jobId']
            logger.info(f"AWS Batch job submitted: {job_id}")
            
        except ImportError:
            logger.warning("boto3 not installed. Install with: pip install boto3")
            logger.info("Generating deployment configuration instead")
            job_id = f"simulated-aws-{pipeline_id}"
        
        result = {
            'job_id': job_id,
            'platform': 'aws',
            'region': self.region,
            'status': 'submitted',
            'job_definition': job_def,
            'monitoring_url': f"https://console.aws.amazon.com/batch/home?region={self.region}#/jobs"
        }
        
        return result
    
    def _deploy_gcp(
        self,
        pipeline_id: int,
        data_location: str,
        output_location: str,
        instance_type: Optional[str],
        docker_image: str,
        **kwargs
    ) -> Dict:
        """Deploy on Google Cloud using Cloud Life Sciences API."""
        logger.info("Deploying on Google Cloud Platform")
        
        # Default instance type
        if not instance_type:
            instance_type = "n1-standard-8"  # 8 vCPUs, 30 GB RAM
        
        # Create job definition
        job_def = {
            'pipeline_id': pipeline_id,
            'platform': 'gcp',
            'region': self.region,
            'machine_type': instance_type,
            'docker_image': docker_image,
            'data_location': data_location,
            'output_location': output_location
        }
        
        # For actual GCP deployment (requires google-cloud libraries)
        try:
            from google.cloud import lifesciences_v2beta
            
            client = lifesciences_v2beta.WorkflowsServiceV2BetaClient()
            
            # Define pipeline
            pipeline = {
                'actions': [
                    {
                        'imageUri': docker_image,
                        'commands': [
                            'python', '/raptor/scripts/01_run_all_pipelines_python.py',
                            '--pipelines', str(pipeline_id),
                            '--data', data_location,
                            '--output', output_location
                        ]
                    }
                ],
                'resources': {
                    'regions': [self.region],
                    'virtualMachine': {
                        'machineType': instance_type,
                        'bootDiskSizeGb': 100
                    }
                }
            }
            
            # Run pipeline
            project = kwargs.get('project_id', os.getenv('GCP_PROJECT'))
            parent = f"projects/{project}/locations/{self.region}"
            
            operation = client.run_pipeline(parent=parent, pipeline=pipeline)
            job_id = operation.name
            
            logger.info(f"GCP pipeline submitted: {job_id}")
            
        except ImportError:
            logger.warning("google-cloud-lifesciences not installed")
            logger.info("Install with: pip install google-cloud-lifesciences")
            logger.info("Generating deployment configuration instead")
            job_id = f"simulated-gcp-{pipeline_id}"
        
        result = {
            'job_id': job_id,
            'platform': 'gcp',
            'region': self.region,
            'status': 'submitted',
            'job_definition': job_def,
            'monitoring_url': f"https://console.cloud.google.com/lifesciences/operations"
        }
        
        return result
    
    def _deploy_azure(
        self,
        pipeline_id: int,
        data_location: str,
        output_location: str,
        instance_type: Optional[str],
        docker_image: str,
        **kwargs
    ) -> Dict:
        """Deploy on Azure using Container Instances or Batch."""
        logger.info("Deploying on Microsoft Azure")
        
        # Default instance type
        if not instance_type:
            instance_type = "Standard_D8s_v3"  # 8 vCPUs, 32 GB RAM
        
        # Create job definition
        job_def = {
            'pipeline_id': pipeline_id,
            'platform': 'azure',
            'region': self.region,
            'vm_size': instance_type,
            'docker_image': docker_image,
            'data_location': data_location,
            'output_location': output_location
        }
        
        # For actual Azure deployment (requires azure-mgmt libraries)
        try:
            from azure.mgmt.containerinstance import ContainerInstanceManagementClient
            from azure.identity import DefaultAzureCredential
            
            credential = DefaultAzureCredential()
            subscription_id = kwargs.get('subscription_id', os.getenv('AZURE_SUBSCRIPTION_ID'))
            
            client = ContainerInstanceManagementClient(credential, subscription_id)
            
            # Define container group
            container_group = {
                'location': self.region,
                'containers': [
                    {
                        'name': f'raptor-pipeline-{pipeline_id}',
                        'image': docker_image,
                        'command': [
                            'python', '/raptor/scripts/01_run_all_pipelines_python.py',
                            '--pipelines', str(pipeline_id),
                            '--data', data_location,
                            '--output', output_location
                        ],
                        'resources': {
                            'requests': {
                                'cpu': 8,
                                'memoryInGB': 32
                            }
                        }
                    }
                ],
                'osType': 'Linux',
                'restartPolicy': 'Never'
            }
            
            resource_group = kwargs.get('resource_group', 'raptor-resources')
            container_group_name = f'raptor-pipeline-{pipeline_id}'
            
            # Create container group
            operation = client.container_groups.begin_create_or_update(
                resource_group,
                container_group_name,
                container_group
            )
            
            job_id = f"{resource_group}/{container_group_name}"
            logger.info(f"Azure container created: {job_id}")
            
        except ImportError:
            logger.warning("azure-mgmt-containerinstance not installed")
            logger.info("Install with: pip install azure-mgmt-containerinstance azure-identity")
            logger.info("Generating deployment configuration instead")
            job_id = f"simulated-azure-{pipeline_id}"
        
        result = {
            'job_id': job_id,
            'platform': 'azure',
            'region': self.region,
            'status': 'submitted',
            'job_definition': job_def,
            'monitoring_url': f"https://portal.azure.com/#blade/HubsExtension/BrowseResource/resourceType/Microsoft.ContainerInstance%2FcontainerGroups"
        }
        
        return result
    
    def check_status(self, job_id: str) -> Dict:
        """
        Check status of cloud job.
        
        Parameters
        ----------
        job_id : str
            Job identifier returned by deploy_pipeline
        
        Returns
        -------
        dict
            Job status information
        """
        logger.info(f"Checking status for job: {job_id}")
        
        if self.platform == 'aws':
            return self._check_aws_status(job_id)
        elif self.platform == 'gcp':
            return self._check_gcp_status(job_id)
        elif self.platform == 'azure':
            return self._check_azure_status(job_id)
    
    def _check_aws_status(self, job_id: str) -> Dict:
        """Check AWS Batch job status."""
        try:
            import boto3
            
            batch_client = boto3.client('batch', region_name=self.region)
            response = batch_client.describe_jobs(jobs=[job_id])
            
            if response['jobs']:
                job = response['jobs'][0]
                status = {
                    'job_id': job_id,
                    'status': job['status'],
                    'created_at': str(job.get('createdAt', '')),
                    'started_at': str(job.get('startedAt', '')),
                    'stopped_at': str(job.get('stoppedAt', ''))
                }
            else:
                status = {'job_id': job_id, 'status': 'not_found'}
            
        except Exception as e:
            logger.error(f"Error checking AWS status: {e}")
            status = {'job_id': job_id, 'status': 'error', 'error': str(e)}
        
        return status
    
    def _check_gcp_status(self, job_id: str) -> Dict:
        """Check GCP operation status."""
        try:
            from google.cloud import lifesciences_v2beta
            
            client = lifesciences_v2beta.WorkflowsServiceV2BetaClient()
            operation = client.get_operation(name=job_id)
            
            status = {
                'job_id': job_id,
                'status': 'running' if not operation.done else 'completed',
                'done': operation.done
            }
            
        except Exception as e:
            logger.error(f"Error checking GCP status: {e}")
            status = {'job_id': job_id, 'status': 'error', 'error': str(e)}
        
        return status
    
    def _check_azure_status(self, job_id: str) -> Dict:
        """Check Azure container status."""
        try:
            from azure.mgmt.containerinstance import ContainerInstanceManagementClient
            from azure.identity import DefaultAzureCredential
            
            credential = DefaultAzureCredential()
            subscription_id = os.getenv('AZURE_SUBSCRIPTION_ID')
            
            client = ContainerInstanceManagementClient(credential, subscription_id)
            
            resource_group, container_group_name = job_id.split('/')
            container_group = client.container_groups.get(resource_group, container_group_name)
            
            status = {
                'job_id': job_id,
                'status': container_group.provisioning_state,
                'state': container_group.containers[0].instance_view.current_state.state if container_group.containers else 'unknown'
            }
            
        except Exception as e:
            logger.error(f"Error checking Azure status: {e}")
            status = {'job_id': job_id, 'status': 'error', 'error': str(e)}
        
        return status
    
    def download_results(
        self,
        output_location: str,
        local_dir: str
    ):
        """
        Download results from cloud storage.
        
        Parameters
        ----------
        output_location : str
            Cloud storage location
        local_dir : str
            Local directory to download to
        """
        logger.info(f"Downloading results from {output_location} to {local_dir}")
        
        local_path = Path(local_dir)
        local_path.mkdir(parents=True, exist_ok=True)
        
        if self.platform == 'aws':
            self._download_from_s3(output_location, local_dir)
        elif self.platform == 'gcp':
            self._download_from_gcs(output_location, local_dir)
        elif self.platform == 'azure':
            self._download_from_azure(output_location, local_dir)
    
    def _download_from_s3(self, s3_location: str, local_dir: str):
        """Download from AWS S3."""
        try:
            subprocess.run(
                ['aws', 's3', 'sync', s3_location, local_dir],
                check=True
            )
            logger.info(f"Downloaded from S3: {s3_location}")
        except Exception as e:
            logger.error(f"Error downloading from S3: {e}")
    
    def _download_from_gcs(self, gcs_location: str, local_dir: str):
        """Download from Google Cloud Storage."""
        try:
            subprocess.run(
                ['gsutil', '-m', 'rsync', '-r', gcs_location, local_dir],
                check=True
            )
            logger.info(f"Downloaded from GCS: {gcs_location}")
        except Exception as e:
            logger.error(f"Error downloading from GCS: {e}")
    
    def _download_from_azure(self, azure_location: str, local_dir: str):
        """Download from Azure Blob Storage."""
        try:
            subprocess.run(
                ['az', 'storage', 'blob', 'download-batch',
                 '--source', azure_location,
                 '--destination', local_dir],
                check=True
            )
            logger.info(f"Downloaded from Azure: {azure_location}")
        except Exception as e:
            logger.error(f"Error downloading from Azure: {e}")


if __name__ == '__main__':
    print("RAPTOR Cloud Integration")
    print("=======================")
    print("\nDeploy RAPTOR pipelines on cloud platforms.")
    print("\nSupported platforms:")
    print("  • AWS (Amazon Web Services)")
    print("  • GCP (Google Cloud Platform)")
    print("  • Azure (Microsoft Azure)")
    print("\nUsage:")
    print("  from raptor.cloud_integration import CloudRunner")
    print("  runner = CloudRunner(platform='aws')")
    print("  runner.deploy_pipeline(pipeline_id=3, data_location='s3://...')")
