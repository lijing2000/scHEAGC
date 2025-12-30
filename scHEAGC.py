import opt
from encoder import *


class scHEAGC(nn.Module):
    def __init__(self, vqvae1, vqvae2, gae1, gae2, n_node=None):
        super(scMIC, self).__init__()

        self.vqvae1 = vqvae1
        self.vqvae2 = vqvae2

        self.gae1 = gae1
        self.gae2 = gae2

        self.a1 = Parameter(nn.init.constant_(torch.zeros(n_node, opt.args.n_z), 0.5), requires_grad=True)  # Z_ae, Z_igae
        self.a2 = Parameter(nn.init.constant_(torch.zeros(n_node, opt.args.n_z), 0.8), requires_grad=True)  # Z_ae, Z_igae
        self.alpha = Parameter(torch.zeros(1))   # ZG, ZL

        self.cluster_centers1 = Parameter(torch.Tensor(opt.args.n_clusters, opt.args.n_z), requires_grad=True)
        self.cluster_centers2 = Parameter(torch.Tensor(opt.args.n_clusters, opt.args.n_z), requires_grad=True)
        torch.nn.init.xavier_normal_(self.cluster_centers1.data)
        torch.nn.init.xavier_normal_(self.cluster_centers2.data)
        self.q_distribution1 = q_distribution(self.cluster_centers1)
        self.q_distribution2 = q_distribution(self.cluster_centers2)

        self.label_contrastive_module = nn.Sequential(
            nn.Linear(n_node, opt.args.n_clusters),
            nn.Softmax(dim=1)
        )

    
    def emb_fusion(self, adj, z_qua, z_igae):
        z_i1 = self.a1 * z_qua + (1 - self.a1) * z_igae
        z_l1 = torch.spmm(adj, z_i1)
    
        z_i2 = self.a2 * z_l1 + (1 - self.a2) * z_qua
        z_l2 = torch.spmm(adj, z_i2)
    
        s = torch.mm(z_l2, z_l2.t())
        s = F.softmax(s, dim=1)
        z_g = torch.mm(s, z_l2)
        z_tilde = self.alpha * z_g + z_l2
        return z_tilde

    def forward(self, x1, adj1, x2, adj2, pretrain=False):

        z_qua1,_ = self.vqvae1.encoder(x1)
        z_qua2,_ = self.vqvae2.encoder(x2)

        # node embedding encoded by IGAE
        z_igae1, a_igae1 = self.gae1.encoder(x1, adj1)
        z_igae2, a_igae2 = self.gae2.encoder(x2, adj2)

        z1 = self.emb_fusion(adj1, z_qua1, z_igae1)
        z2 = self.emb_fusion(adj2, z_qua2, z_igae2)

        z1_tilde = self.label_contrastive_module(z1.T)
        z2_tilde = self.label_contrastive_module(z2.T)

        cons = [z1, z2, z1_tilde, z2_tilde]


        x_hat1 = self.vqvae1.decoder(z1)
        x_hat2 = self.vqvae2.decoder(z2)


        z_hat1, z_adj_hat1 = self.gae1.decoder(z1, adj1)
        a_hat1 = a_igae1 + z_adj_hat1

        z_hat2, z_adj_hat2 = self.gae2.decoder(z2, adj2)
        a_hat2 = a_igae2 + z_adj_hat2

        if not pretrain:
            # the soft assignment distribution Q
            Q1 = self.q_distribution1(z1, z_qua1, z_igae1)
            Q2 = self.q_distribution2(z2, z_qua2, z_igae2)
        else:
            Q1, Q2 = None, None
           
        return x_hat1, z_hat1, a_hat1, x_hat2, z_hat2, a_hat2, Q1, Q2, z1, z2, cons
