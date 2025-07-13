import streamlit as st
import requests, smtplib
from email.mime.text import MIMEText
from datetime import datetime, timedelta
from typing import List, Dict
from Bio import Entrez
from langchain_openai import AzureChatOpenAI

# 🔐 CONFIG — replace with your own
AZURE_ENDPOINT   = "https://YOUR-ENDPOINT.cognitiveservices.azure.com/"
AZURE_API_KEY    = "YOUR-KEY"
CHAT_DEPLOY      = "gpt-4o"
Entrez.email     = "your@email.com"
SMTP_EMAIL       = "your@gmail.com"
SMTP_PASSWORD    = "your-app-password"  # NOT your Gmail password

llm = AzureChatOpenAI(
    azure_endpoint     = AZURE_ENDPOINT,
    api_key            = AZURE_API_KEY,
    azure_deployment   = CHAT_DEPLOY,
    openai_api_version = "2024-12-01-preview",
    temperature        = 0,
)

def fetch_pubmed(topic: str, max_results=10, days=None):
    params = dict(db="pubmed", term=topic, retmax=max_results)
    if days:
        since = (datetime.utcnow() - timedelta(days=days)).strftime("%Y/%m/%d")
        params.update(dict(mindate=since, datetype="pdat"))
    ids = Entrez.read(Entrez.esearch(**params))["IdList"]
    if not ids:
        return []
    docs = Entrez.read(Entrez.esummary(db="pubmed", id=",".join(ids)))
    return [
        dict(title=doc.get("Title", "")[:250],
             date=doc.get("PubDate", "N/A"),
             url=f"https://pubmed.ncbi.nlm.nih.gov/{doc.get('Id')}/")
        for doc in docs
    ]

def fetch_trials(condition, keyword, max_results=5):
    r = requests.get("https://clinicaltrials.gov/api/v2/studies", params={
        "query.cond": condition,
        "query.term": keyword,
        "pageSize":   max_results,
        "format":     "json"
    }, timeout=20)
    r.raise_for_status()
    trials = r.json().get("studies", [])
    results = []
    for s in trials:
        ps = s["protocolSection"]
        ident = ps["identificationModule"]
        status = ps["statusModule"]
        phase = ps.get("designModule", {}).get("phaseList", {}).get("phase", ["N/A"])[0]
        results.append({
            "title": ident.get("officialTitle") or ident.get("briefTitle"),
            "nct": ident["nctId"],
            "phase": phase,
            "start": status.get("startDateStruct", {}).get("startDate", "N/A"),
            "url": f"https://clinicaltrials.gov/study/{ident['nctId']}"
        })
    return results

def summarize(text: str) -> str:
    prompt = (
        "Summarize the following item for a weekly oncology digest in ≤30 words, "
        "focusing on clinical relevance. Do not fabricate or overstate.\n\n"
        f"Item:\n{text}\n\nBullet:"
    )
    return llm.invoke(prompt).content.strip()

def generate_digest(papers, trials, topic):
    today = datetime.now().strftime("%d %b %Y")
    lines = [
        f"Subject: Weekly Oncology Digest – {topic} ({today})",
        "",
        "Dear Dr. Bukavina,",
        "",
        f"Here is your briefing on *{topic}*:",
        "",
        "📚 Latest Papers:"
    ]
    if papers:
        for p in papers:
            bullet = summarize(p['title'])
            lines.append(f"• {bullet}")
            lines.append(f"  ({p['date']}) {p['url']}")
    else:
        lines.append("  – No new articles found.")

    lines.append("")
    lines.append("🧪 Clinical Trials:")
    if trials:
        for t in trials:
            bullet = summarize(t['title'])
            lines.append(f"• {bullet} [{t['phase']}]")
            lines.append(f"  NCT {t['nct']} — Start {t['start']}")
            lines.append(f"  {t['url']}")
    else:
        lines.append("  – No matching trials found.")

    lines += ["", "Warm regards,", "Your AI Assistant"]
    return "\n".join(lines)

def send_email(to_email, subject, body):
    msg = MIMEText(body)
    msg["Subject"] = subject
    msg["From"] = SMTP_EMAIL
    msg["To"] = to_email

    with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
        server.login(SMTP_EMAIL, SMTP_PASSWORD)
        server.send_message(msg)

# ── STREAMLIT UI ───────────────────────────────────────────────
st.title("🧬 Weekly Oncology Digest Generator")
topic = st.text_input("Topic (e.g. ctDNA and bladder cancer)", "ctDNA and bladder cancer")
email = st.text_input("Recipient Email")
max_papers = st.slider("Max Papers", 1, 50, 8)
max_trials = st.slider("Max Trials", 1, 20, 5)
all_time = st.checkbox("Include all-time PubMed results?", value=False)

if st.button("Generate & Send"):
    with st.spinner("Fetching research..."):
        condition = "bladder cancer"
        keyword = "ctDNA" if "ctDNA" in topic else topic
        papers = fetch_pubmed(topic, max_papers, None if all_time else 7)
        trials = fetch_trials(condition, keyword, max_trials)
        digest = generate_digest(papers, trials, topic)

    st.success("Digest generated!")
    st.text_area("Preview", digest, height=400)

    if email:
        send_email(email, f"Weekly Oncology Digest – {topic}", digest)
        st.success(f"Email sent to {email}")
    else:
        st.warning("Enter an email to send.")
